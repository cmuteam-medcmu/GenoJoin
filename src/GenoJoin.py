import os
import psutil
import threading
import numpy as np
from time import time
from extract_regions import Extract_region
from sentence2vec import WordDatabase

def track_usage(stop_flag, results):
    p = psutil.Process(os.getpid())
    max_cpu = max_mem = 0
    while not stop_flag.is_set():
        cpu = p.cpu_percent(interval=1)
        mem = p.memory_info().rss / (1024**2)
        max_cpu = max(max_cpu, cpu)
        max_mem = max(max_mem, mem)
    results["cpu"] = max_cpu
    results["mem"] = max_mem

def monitor_function(func, *args, **kwargs):
    stop_flag = threading.Event()
    results = {}
    t = threading.Thread(target=track_usage, args=(stop_flag, results))
    t.start()

    # run target function
    func(*args, **kwargs)

    # stop monitoring
    stop_flag.set()
    t.join()
    print(f"Max CPU: {(np.round(results['cpu'])/100):.2f} | Max Mem: {results['mem']:.2f} MB")


def main():
    dirpath = '/mnt/hdd/public/share_resource/CRU_cfDNA/Code_session/McFeatures/Term_project/dv_vcf'

    # Extract 
    all_var, uniq_var = Extract_region(dirpath, "chr12", 1, 1_000_000)
    print(f'Total variants : {len(all_var)}')
    print(f'Total unique   : {len(uniq_var)}')

    w = WordDatabase(uniq_var)
    w.CreateDB(10)

    first_items = []
    for sublist in all_var:
        first_items.append(sublist[0])

    w.QueryDB(first_items)


start_time = time()
monitor_function(main)
print(f'\n>> Run Time: {(time() - start_time):.2f} sec')