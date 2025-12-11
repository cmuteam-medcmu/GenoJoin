import logging
import os
from datetime import datetime
from logging import FileHandler, StreamHandler
from time import time

import pandas as pd


def CheckExist(path: str, key: str) -> bool:
    if key == "dir":
        return os.path.isdir(path)
    elif key == "file":
        return os.path.isfile(path)


def TimeStamp(key: str, start_time: time) -> str:
    stamp_headline = {
        "join": "Join Variants",
        "total": "Total Time",
    }
    end_time = time() - start_time

    return f"{stamp_headline[key]}: {end_time:.2f} sec"


class ThreadsManager:
    def __init__(self, cfg):
        self.inputdir = cfg.raw["input_vcf"]
        self.outname = cfg.raw["output_name"]
        self.dbdir = cfg.raw["outdir"]
        self.qual = cfg.raw["min_qual"]
        self.depth = cfg.raw["min_depth"]
        self.gq = cfg.raw["min_gq"]
        self.threads = cfg.raw["threads"]
        self.loc = pd.read_csv(
            cfg.raw["regions"], header=None, sep="\t", names=["CHROM", "START", "END"]
        )

    def singleThreads(self):
        Chrom_List = []
        for pos in range(self.loc.shape[0]):
            Chrom_List.append(
                [
                    [
                        self.inputdir,
                        self.outname,
                        self.dbdir,
                        self.loc["CHROM"][pos],
                        int(self.loc["START"][pos]),
                        int(self.loc["END"][pos]),
                        self.qual,
                        self.depth,
                        self.gq,
                    ]
                ]
            )
        return 1, Chrom_List

    def multiThreads(self):
        lines = len(self.loc)
        Wind_List = []
        sep = lines // self.threads
        integer = lines % self.threads
        group = []
        for i in range(self.threads):
            if integer > 0:
                element = sep + 1
            elif integer <= 0:
                element = sep
            integer -= 1
            group.append(element)

        num_start = 0
        num_end = group[0] - 1

        for i in range(self.threads):
            Chrom_List = []
            for pos in range(self.loc.shape[0]):
                if pos >= num_start and pos <= num_end:
                    Chrom_List.append(
                        [
                            self.inputdir,
                            self.outname,
                            self.dbdir,
                            self.loc["CHROM"][pos],
                            int(self.loc["START"][pos]),
                            int(self.loc["END"][pos]),
                            self.qual,
                            self.depth,
                            self.gq,
                        ]
                    )
                else:
                    continue
            num_start = num_end + 1
            if i < self.threads - 1:
                num_end = num_end + group[i + 1]
            else:
                num_end = lines
            Wind_List.append(Chrom_List)
        return len(Wind_List), Wind_List

    def main(self):
        if self.threads > 1:
            dataList = self.multiThreads()
        else:
            dataList = self.singleThreads()

        return dataList


class Logger:
    def __init__(self, outdir):
        now = datetime.now().strftime("%Y-%m-%d")
        self.log_path = os.path.join(outdir, f"monitor-{now}.log")

        formatter = logging.Formatter("%(asctime)s | %(message)s", "%Y-%m-%d %H:%M:%S")

        self.logger = logging.getLogger(f"GenoJoin")

        # File handler
        fh = FileHandler(self.log_path)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)

        # Console handler
        ch = StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)

        # Register handlers
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)

        self.logo: str = """================================================
  _____                       _       _       
 / ____|                     | |     (_)      
| |  __  ___ _ __   ___      | | ___  _ _ __  
| | |_ |/ _ \ '_ \ / _ \ _   | |/ _ \| | '_ \ 
| |__| |  __/ | | | (_) | |__| | (_) | | | | |
 \_____|\___|_| |_|\___/ \____/ \___/|_|_| |_|
                                            
================================================
          Songphon Sutthitthasakul  

        CMUTEAM, Faculty of Medicine, 
            Chiang Mai University
================================================\n"""

    def run(self):

        print(self.logo)

        with open(self.log_path, "a") as log:
            log.write(f"{self.logo}\n")
            log.close()

        return self.logger


def SummaryVariants(all_vars: list, uniq_vars: list, samples: list, filter: int) -> str:
    return f"Total samples     : {len(samples)}\n                    | Total variants    : {len(all_vars)}\n                    | Total unique      : {len(uniq_vars)}\n                    | Filtered variants : {filter}"


def FormatTransform(format: str, variant: str) -> list:
    form = format.split(":")
    var = variant.split(":")
    new_var = []

    if len(form) < 4:
        return new_var, False

    new_var.append(var[form.index("GT")])

    if "DP" in format:
        try:
            DP = form.index("DP")
        except ValueError:
            DP = -1
    elif "MIN_DP" in format:
        try:
            DP = form.index("MIN_DP")
        except ValueError:
            DP = -1
    else:
        DP = -1

    if DP == -1:
        return new_var, False

    new_var.append(var[form.index("GQ")])
    new_var.append(var[DP])

    try:
        AD = form.index("AD")
    except ValueError:
        AD = -1

    if AD == -1:
        new_var.append(f"{var[DP]},0")
    else:
        new_var.append(var[form.index("AD")])

    new_var.append(var[form.index("PL")])

    return new_var, True


def get_prefix(s: str) -> str:
    parts = s.split(" ")
    parts[3] = parts[3][0]
    return " ".join(parts[:3])


def get_pattern(s: str) -> str:
    parts = s.split(":")
    return parts[0]


def RemoveRef(A: list, B: list) -> list:
    prefix_A = [get_prefix(x[0]) for x in A]
    filtered_B = [x for x in B if get_prefix(x[0]) not in prefix_A]
    return filtered_B


def RemoveNoise(A: list, B: list) -> list:
    prefix_A = [get_prefix(x) for x in A]
    filtered_B = [x for x in B if get_prefix(x[0]) in prefix_A]
    return filtered_B
