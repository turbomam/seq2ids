# import pprint
# from gffutils.interface import FeatureDB
# from gffutils.interface import Feature
import ast
import datetime
import sqlite3
from os import listdir
from os.path import isfile, join, splitext
from typing import List

import gffutils
import pandas as pd

pd.set_option("display.max_columns", None)

# todo iterate over blast hit coordinates
# todo determine allen interval type
# todo special deletion flank handling

gffcols = {
    "id": 0,
    "seqid": 1,
    "source": 2,
    "featuretype": 3,
    "start": 4,
    "end": 5,
    "score": 6,
    "strand": 7,
    "frame": 8,
    "attributes": 9,
    "extra": 10,
    "bin": 11,
}

ctf_fp = "ranges_to_features.tsv"

skippers = ["attributes", "extra"]

gff_dir = "../gff/"

# ---

seq2ids_db_fp = "../target/seq2ids.db"

# ---

all_targets = list(gffcols.keys())

targets = [x for x in all_targets if x not in skippers]


class Features:
    def __init__(self):
        self.something = None
        self.lod: List = []

    def do_something(self, something):
        self.something = something

    def feature_tuple_to_dict(self, feature_tuple: tuple) -> dict:
        f_keys = list(gffcols.keys())
        as_dict = dict(zip(f_keys, feature_tuple))
        return as_dict

    def gff_to_db(self, gff_file: str) -> gffutils.FeatureDB:
        db = gffutils.create_db(
            gff_file, ":memory:", merge_strategy="create_unique", keep_order=True, force=True
        )
        # db = gffutils.FeatureDB('myGFF.db')
        db.analyze()
        return db

    def return_frame(self):
        frame = pd.DataFrame(self.lod)
        return frame

    def coords_to_frame(
            self, db: gffutils.FeatureDB, start: int, end: int
    ) -> pd.DataFrame:

        extract = db.region(start=start, end=end)
        for i in extract:
            it = i.astuple()
            itd = self.feature_tuple_to_dict(it)
            it_ft = itd["featuretype"]
            if it_ft != "region":
                rearranged = {key: itd[key] for key in itd.keys() if key in targets}
                attributes = ast.literal_eval(itd["attributes"])

                # assume a single dict of single-element lists?
                flattened = {key: value[0] for key, value in attributes.items()}

                merged = rearranged.copy()
                for key, value in flattened.items():
                    merged[key] = value

                self.lod.append(merged)
        # db.close()


fancy_instance = Features()
seq2ids_db_conn = sqlite3.connect(seq2ids_db_fp)
seq2ids_db_curs = seq2ids_db_conn.cursor()

gff_file_list = [f for f in listdir(gff_dir) if isfile(join(gff_dir, f))]
gff_file_list.sort()

df_list = []

for gff_file in gff_file_list:
    genome_name = splitext(gff_file)[0]
    print(genome_name)
    ts = datetime.datetime.now().timestamp()
    readable = datetime.datetime.fromtimestamp(ts).isoformat()
    print(readable)
    blast_hits_q = f"SELECT * FROM blast_results where sacc = '{genome_name}'"
    blast_hits_res = pd.read_sql_query(blast_hits_q, seq2ids_db_conn)
    bhr_lod = blast_hits_res.to_dict(orient="records")
    gff_fp = join(gff_dir, gff_file)
    gff_db = fancy_instance.gff_to_db(gff_fp)
    for bh in bhr_lod:
        fancy_instance.coords_to_frame(gff_db, bh["sstart"], bh["send"])
        from_dump = fancy_instance.return_frame()
        from_dump.insert(0, "qacc", bh["qacc"])
        from_dump.insert(1, "qstart", bh["qstart"])
        from_dump.insert(2, "qend", bh["qend"])
        # there are id and ID columns. I have not seen a case in which their content has been different.
        # sqlite complains about columns with the same names case-insensitive
        from_dump.drop(columns=['ID'], inplace=True)
        df_list.append(from_dump.columns)
#         from_dump.dropna(axis=1, how='all', inplace=True)
#         df_list.append(from_dump)
#
# concat = pd.concat(df_list)
# concat.to_csv("concat.tsv", sep="\t", index=False)
# concat.to_sql(name="features", con=seq2ids_db_conn, index=False)

# todo add to database... be prepared for different columns... aggregate as one frame first?

# # # evalue | bitscore | "length" | pident | qacc | qcovhsp | qcovs | qstart | qend | qseqid | sacc | sallacc | sallgi |
# # # sallseqid | salltitles | sblastnames | scomnames | sstart | send | sgi | sseqid | staxids | stitle | sskingdoms |
# # # ssciname | gapopen | mismatch

flat_list = [item for sublist in df_list for item in sublist]

print(flat_list)