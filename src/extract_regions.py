import os
import tabix
from time import time
from Bio import SeqIO


def Extract_region(dirpath, chr, start, end):
    start_time = time()
    fasta_path = '/mnt/sas/ref/hg38/v0/Homo_sapiens_assembly38.fasta'
    variants = []
    oo = 0
    for path in os.listdir(dirpath):
        if path.endswith("g.vcf.gz"):
            tbx = tabix.open(f"{dirpath}/{path}")

            records = tbx.query(chr, start, end)
            for vix, var in enumerate(records):
                
                if var[4] == '<*>' or var[4] == '':
                    variants.append([f'{var[0]} {var[1]} {var[3]} {var[3]}', oo])
                    
                    continue
                elif ',<*>' in var[4]:
                    var[4] = var[4][:-4]
                
                pos = int(var[1])

                if var[4].count(',') > 0:
                    for allele in var[4].split(','):
                        if len(var[3]) == len(allele) and len(var[3]) > 1:
                            for zix, (zet1, zet2) in enumerate(zip(var[3],allele)):
                                var[1] = pos + zix
                                variants.append([f'{var[0]} {var[1]} {zet1} {zet2}', oo]) 
                                print([f'{var[0]} {var[1]} {zet1} {zet2}', oo])
                        else:
                            variants.append([f'{var[0]} {var[1]} {var[3]} {allele}', oo]) 

                else:
                    if len(var[3]) == len(var[4]) and len(var[3]) > 1:
                        for zix, (zet1, zet2) in enumerate(zip(var[3],var[4])):
                            var[1] = pos + zix
                            variants.append([f'{var[0]} {var[1]} {zet1} {zet2}', oo])
                    else:
                        variants.append([f'{var[0]} {var[1]} {var[3]} {var[4]}', oo])

            oo += 1

    end_time = time() - start_time
    print(f'>> Extract Regions: {end_time:.2f} sec')
    return variants, list(dict.fromkeys(item[0] for item in variants))