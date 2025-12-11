import os
from time import time

import tabix

from utils import FormatTransform, TimeStamp


def Extract_region(dirpath, chr, start, end, qual, depth, gq):
    start_time = time()
    germline = []
    variants = []
    for path in sorted(os.listdir(dirpath)):
        if path.endswith("norm.vcf.gz"):

            tbx = tabix.open(f"{dirpath}/{path}")

            records = tbx.query(chr, start, end)
            for vix, var in enumerate(records):

                # Format must be (GT:GQ:DP:AD:AF:PL)
                eachVar = FormatTransform(var[8], var[9])
                eachVar.append(var[5])

                if int(eachVar[1]) < gq:
                    continue

                ADs = list(map(int, eachVar[3].split(",")))
                if var[4] == "<*>" or var[4] == "<NON_REF>" or var[4] == "":
                    var[4] = var[3]
                    if int(eachVar[2]) < depth:
                        eachVar[0] = "./."
                    else:
                        eachVar[0] = "0/0"

                    germline.append(
                        [
                            f"{var[0]} {var[1]} {var[3]} {var[4]}",
                            ":".join(eachVar),
                            path.split("_")[0],
                        ]
                    )

                else:
                    if float(var[5]) >= qual:
                        AD = ADs[1] / (ADs[0] + ADs[1])
                        if AD > 0.2 and AD <= 0.8:
                            eachVar[0] = "0/1"
                        elif AD <= 0.2:
                            continue
                        elif AD > 0.8:
                            eachVar[0] = "1/1"

                        variants.append(
                            [
                                f"{var[0]} {var[1]} {var[3]} {var[4]}",
                                ":".join(eachVar),
                                path.split("_")[0],
                            ]
                        )

    TimeStamp("extract", start_time)
    return (
        germline,
        variants,
        list(dict.fromkeys(item[0] for item in variants)),
        list(dict.fromkeys(item[2] for item in variants)),
    )
