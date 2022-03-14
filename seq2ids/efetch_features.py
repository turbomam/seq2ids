import sqlite3

import numpy as np
import pandas as pd
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "MAM@lbl.gov"
# todo tool?
# todo persistent session?
seq2ids_db_fp = "../target/seq2ids.db"

row_list = []

seq2ids_db_conn = sqlite3.connect(seq2ids_db_fp)
ranges_to_download_q = "select * from ranges_to_download order by sacc"
ranges_to_download_res = pd.read_sql_query(ranges_to_download_q, seq2ids_db_conn)

rtd_lod = ranges_to_download_res.to_dict(orient="records")
for genome_range in rtd_lod:
    print(genome_range)

    handle = Entrez.efetch(
        db="nucleotide",
        id=genome_range["sacc"],
        rettype="gb",
        retmode="text",
        seq_start=genome_range["minpos"],
        seq_stop=genome_range["maxpos"],
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()

    print(f"record.id: {record.id}")

    feature_row = {}
    for i in record.features:
        feature_row["record_id"] = record.id
        feature_row["record_name"] = record.name
        feature_row["record_description"] = record.description
        feature_row["record_dbxrefs"] = '|'.join(record.dbxrefs)

        for rak, rav in record.annotations.items():
            # record_annotations_structured_comment
            # OrderedDict([('Genome-Assembly-Data', OrderedDict([('Assembly Method', 'Velevt v. 1.1.02'), ('Genome Coverage', '80x'), ('Sequencing Technology', 'llumina Solexa')]))])
            if rak not in ["references", "comment", "structured_comment"]:
                if isinstance(rav, list):
                    rav = '|'.join(rav)
                feature_row[f"record_annotations_{rak}"] = rav

        feature_row["feature_id"] = i.id
        feature_row["feature_type"] = i.type
        feature_row["feature_strand"] = i.strand
        # class Bio.SeqFeature.FeatureLocation(start, end, strand=None, ref=None, ref_db=None)
        # feature_row["feature_location"] = str(i.location)
        feature_row["feature_location_start"] = i.location.start
        feature_row["feature_location_end"] = i.location.end
        feature_row["feature_location_strand"] = i.location.strand
        feature_row["feature_location_ref"] = i.location.ref
        feature_row["feature_location_ref_db"] = i.location.ref_db

        feature_row["feature_ref"] = i.ref
        feature_row["feature_ref_db"] = i.ref_db
        feature_row["feature_location_operator"] = i.location_operator

        q = i.qualifiers
        for fqk, fqv in q.items():
            if fqk not in ["references", "comment", "translation"]:
                if isinstance(fqv, list):
                    fqv = '|'.join(fqv)
                feature_row[f"feature_qualifier_{fqk}"] = fqv

        row_list.append(feature_row)

feature_frame = pd.DataFrame(row_list)

feature_frame.replace(r'^\s*$', np.nan, regex=True, inplace=True)

feature_frame.dropna(axis=1, how='all', inplace=True)

feature_frame.to_csv("feature_frame.tsv", sep="\t", index=False)

feature_frame.to_sql("features", seq2ids_db_conn, if_exists="replace")
