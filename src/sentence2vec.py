import os
import subprocess
from collections import Counter
from pathlib import Path

import faiss
import numpy as np
import pandas as pd
from gensim.models import Word2Vec
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from schema import Base, Sample, Variant
from utils import RemoveNoise, RemoveRef


class WordDatabase:
    def __init__(self, sentences_train, region, dbdir):

        # SQLite database
        os.mkdir(f"{dbdir}/GenoJion_DB")

        if os.path.exists(f"{dbdir}/GenoJion_DB/GJ_{region}.db"):
            print(f"Already has GJ_{region}.db")
            os._exit(1)

        self.dbdir = dbdir

        engine = create_engine(
            f"sqlite:///{dbdir}/GenoJion_DB/GJ_{region}.db", echo=False
        )

        # Create table
        Base.metadata.create_all(engine)

        # Session for transactions
        Session = sessionmaker(bind=engine)
        self.session = Session()
        self.sentences_train = sentences_train

        BASE_DIR = Path(__file__).resolve().parent.parents[0]
        self.vcf_header = BASE_DIR / "data" / "vcf" / "vcf_header.txt"

    def createDict(self):
        samples = [
            v[0]
            for v in self.session.query(Sample.name)
            .distinct()
            .order_by(Sample.name)
            .all()
        ]

        samples_dict = {
            "#CHROM": None,
            "POS": None,
            "ID": ".",
            "REF": None,
            "ALT": None,
            "QUAL": 0,
            "FILTER": ".",
            "INFO": ".",
            "FORMAT": "GT:GQ:DP:AD:PL",
        }

        for sample in samples:
            samples_dict[sample] = None

        return samples_dict

    @classmethod
    def sentence_vector(cls, tokens, model):
        vecs = [model.wv[word] for word in tokens if word in model.wv]
        if len(vecs) == 0:
            return np.zeros(model.vector_size)
        return np.mean(vecs, axis=0)

    def compare_ref(self, var: str, ref: str) -> bool:

        parts_var = var.split(" ")
        parts_ref = ref.split(" ")

        # compare first two fields
        return (
            parts_var[0] == parts_ref[0]
            and parts_var[1] == parts_ref[1]
            and parts_var[2][0] == parts_ref[2]
        )

    def CreateDB(self, threads=1):
        # Tokenize
        tokenized_sentences = [s.split() for s in self.sentences_train]

        # Train Word2Vec
        self.model = Word2Vec(
            tokenized_sentences,
            vector_size=300,
            window=5,
            min_count=1,
            workers=threads,
            seed=42,
        )

        sentence_vectors = np.array(
            [self.sentence_vector(t, self.model) for t in tokenized_sentences],
            dtype="float32",
        )

        d = self.model.vector_size  # vector dimension

        faiss.normalize_L2(sentence_vectors)
        self.index = faiss.IndexFlatIP(d)
        self.index.add(sentence_vectors)

        self.session.execute(
            Variant.__table__.insert(),
            [
                {
                    "chrom": v.split(" ")[0],
                    "pos": int(v.split(" ")[1]),
                    "ref": v.split(" ")[2],
                    "alt": v.split(" ")[3],
                }
                for v in self.sentences_train
            ],
        )
        self.session.commit()

    def MapDB(self, sample, sentences_var, sentences_ref):

        # Variants Mapping
        sentences_var_query = [item[0] for item in sentences_var]
        tokenized_var_query = [s.split() for s in sentences_var_query]
        vectors_var_query = np.array(
            [self.sentence_vector(t, self.model) for t in tokenized_var_query],
            dtype="float32",
        )
        faiss.normalize_L2(vectors_var_query)

        distances, indices = self.index.search(vectors_var_query, 1)
        map_data = []
        for row, dataset_idx in enumerate(indices[:, 0]):
            if self.sentences_train[dataset_idx] != sentences_var_query[row]:
                continue
            else:
                map_data.append(
                    {
                        "name": sample,
                        "vid": int(dataset_idx) + 1,
                        "pat": sentences_var[row][1],
                    }
                )

        # References Mapping
        sentences_ref = RemoveRef(sentences_var, sentences_ref)
        sentences_ref = RemoveNoise(self.sentences_train, sentences_ref)
        sentences_ref_query = [item[0] for item in sentences_ref]
        if len(sentences_ref_query) != 0:
            tokenized_ref_query = [s.split() for s in sentences_ref_query]
            vectors_ref_query = np.array(
                [self.sentence_vector(t, self.model) for t in tokenized_ref_query],
                dtype="float32",
            )
            faiss.normalize_L2(vectors_ref_query)

            distances, indices = self.index.search(vectors_ref_query, 1)
            for row, dataset_idx in enumerate(indices[:, 0]):
                if not self.compare_ref(
                    self.sentences_train[dataset_idx], sentences_ref_query[row]
                ):
                    continue
                else:
                    map_data.append(
                        {
                            "name": sample,
                            "vid": int(dataset_idx) + 1,
                            "pat": sentences_ref[row][1],
                        }
                    )

        sampleMap = [Sample(**m) for m in map_data]
        self.session.bulk_save_objects(sampleMap)
        self.session.commit()

    def QueryDB(self, fname, regions):
        ref_dict = self.createDict()

        no_var = "./.:.:.:.,.:.,.,."

        var_ids = (
            self.session.query(Variant)
            .order_by(Variant.pos, Variant.ref, Variant.alt)
            .all()
        )

        vars_list = [v.__dict__ for v in var_ids]
        for v in vars_list:
            v.pop("_sa_instance_state", None)

        join_vars = []

        for vix, var_id in enumerate([v["id"] for v in vars_list]):
            line = ref_dict.copy()

            samples = (
                self.session.query(Sample.name, Sample.pat)
                .filter(Sample.vid == var_id)
                .all()
            )

            samples_list = [{"name": s[0], "pat": s[1]} for s in samples]

            geno_counts = Counter(s["pat"].split(":")[0] for s in samples_list)
            if geno_counts.get("0/1", 0) + geno_counts.get("1/1", 0) <= 1:
                continue
            else:
                line["#CHROM"] = vars_list[vix]["chrom"]
                line["POS"] = vars_list[vix]["pos"]
                line["ID"] = (
                    f'{vars_list[vix]["chrom"]}_{vars_list[vix]["pos"]}_{vars_list[vix]["ref"]}_{vars_list[vix]["alt"]}'
                )
                line["REF"] = vars_list[vix]["ref"]
                line["ALT"] = vars_list[vix]["alt"]
                qual = []
                for samp, pat in ((s["name"], s["pat"]) for s in samples_list):
                    line[samp] = ":".join(pat.split(":")[:-2])
                    qual.append(round(float(pat.split(":")[-1]), 2))

                line["QUAL"] = np.median(qual)
                join_vars.append(line)

        self.filter_var = len(join_vars)
        df = pd.DataFrame(join_vars)
        df.iloc[:, 9:] = df.iloc[:, 9:].fillna(no_var)
        df.to_csv(
            f"{self.dbdir}/{fname}_{regions}.vcf", sep="\t", header=True, index=False
        )

        try:
            with open(f"{self.dbdir}/{fname}_{regions}_h.vcf", "w") as out_file:
                subprocess.run(
                    ["cat", self.vcf_header, f"{self.dbdir}/{fname}_{regions}.vcf"],
                    stdout=out_file,
                    check=True,
                )

            os.remove(f"{self.dbdir}/{fname}_{regions}.vcf")
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e}")
