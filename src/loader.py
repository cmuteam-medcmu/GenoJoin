import argparse
import os
from dataclasses import dataclass, field

import yaml
from utils import CheckExist


@dataclass
class Config:
    # Default values
    input_vcf: str = None
    min_qual: int = 30
    min_gq: int = 20
    min_depth: int = 10
    output_name: str = "Geno.out"
    outdir: str = "Geno.DB"
    regions: str = None
    threads: int = 1

    # Store raw dict for debugging if needed
    raw: dict = field(default_factory=dict)


class ArgConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description="GenoJoin: fast joint-variant calling for large populations"
        )

        # add normal CLI args
        self.parser.add_argument("-c", "--config", type=str, help="YAML config file")
        self.parser.add_argument("-i", "--input-vcf", type=str)
        self.parser.add_argument("-q", "--min-qual", type=int)
        self.parser.add_argument("-gq", "--min-gq", type=int)
        self.parser.add_argument("-dp", "--min-depth", type=int)
        self.parser.add_argument("-n", "--output-name", type=str)
        self.parser.add_argument("-o", "--outdir", type=str)
        self.parser.add_argument("-r", "--regions", type=str)
        self.parser.add_argument("-t", "--threads", type=int)

    def load(self):
        # 1) parse initial args
        args = self.parser.parse_args()

        config_dict = {}

        # 2) load YAML if provided
        if args.config:
            with open(args.config, "r") as f:
                yaml_config = yaml.safe_load(f) or {}
            config_dict.update(yaml_config)

        # 3) override with CLI args (only if not None)
        cli_dict = {k: v for k, v in vars(args).items() if v is not None}
        config_dict.update(cli_dict)

        # 4) build final config object
        cfg = Config(**{k: v for k, v in config_dict.items() if hasattr(Config, k)})

        cfg.raw = config_dict  # optional

        if not cfg.input_vcf:
            print("ERROR: (input_vcf) Please add input vcf directory!")
            os._exit(1)

        if not CheckExist(cfg.input_vcf, "dir"):
            print(f"ERROR: (input_vcf) {cfg.input_vcf} does NOT exist!")
            os._exit(1)

        if not CheckExist(cfg.outdir, "dir"):
            print(f"ERROR: (outdir) {cfg.outdir} does NOT exist!")
            os._exit(1)

        if cfg.threads < 1:
            print("ERROR: (threads) Thread must more than 0!")
            os._exit(1)

        if not CheckExist(cfg.regions, "file"):
            print(f"ERROR: (regions) {cfg.regions} does NOT exist!")
            os._exit(1)

        return cfg
