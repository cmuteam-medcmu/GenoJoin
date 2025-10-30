import os
import tabix
import numpy as np
import pandas as pd
from time import time
from pyfaidx import Fasta
from cyvcf2 import VCF,Writer

def get_seq(CHR, STR, END) -> str:
    fasta_path = '/mnt/sas/ref/hg38/v0/Homo_sapiens_assembly38.fasta'
    fasta = Fasta(fasta_path)
    return fasta[CHR][STR-1:END-1].seq

def base_conv(seq):
    base_map = {
        "A" : 1,
        "C" : 2,
        "G" : 3,
        "T" : 4,
    }
    lookup = np.full(128, 0, dtype=int)
    for k, v in base_map.items():
        lookup[ord(k)] = v
    
    return lookup[np.array([ord(c) for c in seq], dtype='uint8')]

def uniq_variants(v):
    def genotyping(base, gt, ref=0)  -> str:
        if gt == '.' and ref == 0:
            return '.'
        elif gt == '.' and ref == 1:
            return base[0]
        else:
            return base[int(gt)]

    uniq_v = []
    for i in v:
        genotype = [i[1]]
        if ',' not in i[2]:
            genotype.append(i[2])
            uniq_v.append((i[0],genotyping(genotype,i[3].split('/')[0]),genotyping(genotype,i[3].split('/')[1])))
        else:
            genotype += [gt for gt in i[2].split(',')]

            for j in i[3].split(','):
                uniq_v.append((i[0],genotyping(genotype,j.split('/')[0]),genotyping(genotype,j.split('/')[1])))

    return uniq_v

class DB:
    def __init__(self, chr, start, end): 
        self.chr = chr
        self.seq = get_seq(chr, start, end)
        self.start = start
        self.end = end - start

    def snp_mapping(self, variants, n):
        self.map = np.full((len(variants)+n, self.end), 0, dtype='uint8')
        if n == 0:
            self.map[0,:] = base_conv(self.seq)
        elif n == 1:
            self.map[0,:] = [0 for i in range(self.end)]

        self.index = {f"{i[0]}_{i[1]}_{i[2]}": idx for idx, i in enumerate(variants)}

        for idx, row in enumerate(variants):
            if len(row[2]) == len(row[1]):
                self.map[idx+n,row[0]-self.start:row[0]-self.start+len(row[2])] = base_conv(row[2])
            # elif len(row[2]) < len(row[1]):
            #     self.map[idx+n,row[0]-self.start:row[0]-self.start+len(row[1])] = np.concatenate((base_conv(row[2]) , np.array([0 for i in range(len(row[1])-len(row[2]))])))


start_time = time()
url = '/mnt/hdd/public/share_resource/CRU_cfDNA/Code_session/McFeatures/Term_project/dv_vcf'

loc_chr = "chr12"
loc_start = 111760000
loc_end = 111820000

n = 0

db_all = []
db_idx = {}

for fpath in sorted(os.listdir(url)):
    if fpath.endswith('.vcf.gz'):
        n += 1
        file = f'{url}/{fpath}'
        file_id = fpath.split('.')[0]
        tb = tabix.open(file)
        records = tb.query(loc_chr, loc_start, loc_end)
        data = [(int(i[1]),i[3],i[4][:-3],i[9].split(':')[0]) for i in records]
        uniq_data = uniq_variants(data)
        db = DB(loc_chr, loc_start, loc_end)
        if n == 1:
            db.snp_mapping(uniq_data, 0)
        else:
            db.snp_mapping(uniq_data, 1)
        db_all.append(db.map)
        db_idx[file_id] = db.index
        # df = pd.DataFrame(db.map)
        # df.columns = [i for i in range(loc_start,loc_end)]


db_stack = np.vstack(db_all)
sorted_idx = np.lexsort(db_stack.T[::-1])
arr_sorted = db_stack[sorted_idx]

print(arr_sorted)

groups_sets = {name: set(values) for name, values in db_idx.items()}
all_elements = set().union(*groups_sets.values())
summary = {}
for elem in sorted(all_elements):
    summary[elem] = [name for name, s in groups_sets.items() if elem in s]

output = {}

for key, groups in summary.items():
    a, b = key.split("_", 1)   # split into two parts
    output.setdefault(a, {})
    output[a].setdefault(b, []).extend(groups)

print(output)

run_time = time() - start_time
print(f"Run time: {run_time:.2f} sec")