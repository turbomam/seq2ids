import gffutils
import pprint
from gffutils.interface import FeatureDB
from gffutils.interface import Feature
import ast
import pandas as pd
import sqlite3

# todo iterate over blast hit coordinates
# todo determine allen interval type
# todo special deletion flank handling

gffcols = {"id": 0,
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
           "bin": 11}

all_targets = list(gffcols.keys())

skippers = ['attributes', 'extra']

targets = [x for x in all_targets if x not in skippers]


def feature_tuple_to_dict(feature_tuple: tuple) -> dict:
    f_keys = list(gffcols.keys())
    as_dict = dict(zip(f_keys, feature_tuple))
    return as_dict


def coords_to_frame(gff_file: str, start: int, end: int) -> pd.DataFrame:
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy="create_unique", keep_order=True)
    # db = gffutils.FeatureDB('myGFF.db')
    db.analyze()
    rearranged_list = []

    extract = db.region(start=start, end=end)
    for i in extract:
        it = i.astuple()
        itd = feature_tuple_to_dict(it)
        it_ft = itd['featuretype']
        if it_ft != 'region':
            rearranged = {key: itd[key] for key in itd.keys() if key in targets}
            attributes = ast.literal_eval(itd['attributes'])

            # assume a single dict of single-element lists?
            flattened = {key: value[0] for key, value in attributes.items()}

            merged = rearranged.copy()
            for key, value in flattened.items():
                merged[key] = value

            rearranged_list.append(merged)

    final_dict = pd.DataFrame(rearranged_list)

    return final_dict


my_genome = 'CP034908'
myGFF = f"../local/{my_genome}.gff"
startloc = 58574
endloc = 58580
ctf_fp = 'ranges_to_features.tsv'

connection = sqlite3.connect("../seq2ids.db")
cursor = connection.cursor()
blast_hits_q = 'SELECT * FROM "seq2ids_blastn-nt_out" where sacc = \'CP034908\''

blast_hits_res = pd.read_sql_query(blast_hits_q, connection)

bhr_lod = blast_hits_res.to_dict(orient='records')
first_record = bhr_lod[0]

ctf = coords_to_frame(myGFF, first_record['sstart'], first_record['send'])

ctf.to_csv(ctf_fp, sep='\t', index=False)

# evalue | bitscore | "length" | pident | qacc | qcovhsp | qcovs | qstart | qend | qseqid | sacc | sallacc | sallgi |
# sallseqid | salltitles | sblastnames | scomnames | sstart | send | sgi | sseqid | staxids | stitle | sskingdoms |
# ssciname | gapopen | mismatch
