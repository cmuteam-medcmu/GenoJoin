import threading
from multiprocessing import Pool
from time import time

import psutil
from tqdm import tqdm

from extract_regions import Extract_region
from loader import ArgConfig
from sentence2vec import WordDatabase
from utils import *


def track_usage(pids, stop_flag, result_dict):
    max_cpu = 0
    max_mem = 0

    while not stop_flag.is_set():
        cpu_total = 0
        mem_total = 0

        for pid in pids:
            try:
                p = psutil.Process(pid)
                cpu_total += p.cpu_percent(interval=0.1)
                mem_total += p.memory_info().rss / (1024 * 1024)  # MB
            except psutil.NoSuchProcess:
                continue

        max_cpu = max(max_cpu, cpu_total)
        max_mem = max(max_mem, mem_total)

    result_dict["max_cpu"] = max_cpu
    result_dict["max_mem"] = max_mem


def main(setList):
    n = len(setList)
    with tqdm(total=n, miniters=n // 10, desc="Processing") as pbar:
        for item in setList:
            st = time()

            (
                inputdir,
                outname,
                dbdir,
                chromosome,
                pos_start,
                pos_end,
                qual,
                depth,
                gq,
            ) = item

            # Extract
            all_ref, all_var, uniq_var, samples = Extract_region(
                inputdir, chromosome, pos_start, pos_end, qual, depth, gq
            )

            if len(uniq_var) == 0:
                pbar.update(1)
                continue

            # Create DB
            w = WordDatabase(uniq_var, f"{chromosome}_{pos_start}_{pos_end}", dbdir)
            w.CreateDB(threads=1)

            # Map DB
            for sample_idx, sample in enumerate(samples):
                sentence_var = [var[:2] for var in all_var if var[2] == sample]
                sentence_ref = [var[:2] for var in all_ref if var[2] == sample]
                w.MapDB(sample, sentence_var, sentence_ref)

            # Query DB
            w.QueryDB(outname, f"{chromosome}_{pos_start}_{pos_end}")

            print("\n")
            print(f"[Region] {chromosome}\t{pos_start}\t{pos_end}")
            summary = SummaryVariants(all_var, uniq_var, samples, w.filter_var)
            join_ts = TimeStamp("join", st)

            pbar.update(1)

    return f"[Region] {chromosome}\t{pos_start}\t{pos_end}\n{summary}\n{join_ts}\n"


if __name__ == "__main__":

    start_time = time()
    cfg = ArgConfig().load()

    logger = Logger(cfg.raw["outdir"]).run()
    logger.info(">> Start Processing")
    threads, dataList = ThreadsManager(cfg).main()

    with Pool(processes=threads) as pool:
        pids = [p.pid for p in pool._pool]
        stop_flag = threading.Event()
        result_dict = {}
        monitor_thread = threading.Thread(
            target=track_usage, args=(pids, stop_flag, result_dict)
        )
        monitor_thread.start()

        results = pool.map(main, dataList)

        stop_flag.set()
        monitor_thread.join()
        pool.close()
        pool.join()

        for r in results:
            logger.info(f"{r}")

        print(
            f"[Usage] CPU: {int(result_dict['max_cpu']/100)} cores | MEM: {result_dict['max_mem']:.1f} MB"
        )

        logger.info(
            f"[Usage] CPU: {int(result_dict['max_cpu']/100)} cores | MEM: {result_dict['max_mem']:.1f} MB"
        )

    total_ts = TimeStamp("total", start_time)
    logger.info(f"{total_ts}")
