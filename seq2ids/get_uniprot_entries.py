import io
import random
import sqlite3
from typing import List

import pandas as pd
import requests

pd.set_option('display.max_columns', None)

seq2ids_db_fp = "target/seq2ids.db"

# sacc_table = 'ranges_to_download'
sacc_table = 'one_best_up'

keep_frac = 1

submission_chunk_size = 300

seq2ids_db_conn = sqlite3.connect(seq2ids_db_fp)

ranges_to_download_q = f"select distinct sacc from {sacc_table} order by sacc"

ranges_to_download_res = pd.read_sql_query(ranges_to_download_q, seq2ids_db_conn)

saccs = ranges_to_download_res['sacc'].to_list()
saccs_len = len(saccs)
keep_num = saccs_len * keep_frac
random_saccs = random.sample(saccs, int(keep_num))


def get_uniprot_sequences(uniprot_ids: List, uniprot_cols: List) -> pd.DataFrame:
    """
    https://www.biostars.org/p/94422/

    Retrieve uniprot sequences based on a list of uniprot sequence identifier.

    For large lists it is recommended to perform batch retrieval.

    documentation which columns are available:
    https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames

    this python script is based on
    https://www.biostars.org/p/67822/

    Parameters:
        uniprot_ids: List, list of uniprot identifier
        uniprot_cols: List, list of display columns (https://www.uniprot.org/help/uniprotkb_column_names and https://www.uniprot.org/docs/dbxref)

    Returns:
        pd.DataFrame, pandas dataframe with uniprot id column and sequence
    """

    joined_cols = ','.join(uniprot_cols)

    url = 'https://www.uniprot.org/uploadlists/'  # This is the webserver to retrieve the Uniprot data
    params = {
        'from': "ACC",
        'to': 'ACC',
        'format': 'tab',
        'query': " ".join(uniprot_ids),
        'columns': joined_cols
    }

    response = requests.get(url, params=params)

    # print(response.request.url)

    response_frame = pd.read_csv(io.BytesIO(response.content), sep='\t')
    response_frame = response_frame.iloc[:, :-1]

    return response_frame


include_fields = ['id', 'entry name', 'genes', 'genes(PREFERRED)', 'genes(ALTERNATIVE)', 'genes(OLN)', 'genes(ORF)',
                  'organism-id', 'organism', 'lineage(ALL)', 'protein names', 'proteome', 'database(BioCyc)',
                  'database(BRENDA)', 'database(CDD)', 'comment(ALTERNATIVE PRODUCTS)',
                  'comment(ERRONEOUS GENE MODEL PREDICTION)', 'comment(ERRONEOUS INITIATION)',
                  'comment(SEQUENCE CAUTION)', 'comment(ENZYME REGULATION)', 'comment(FUNCTION)',
                  'ec', 'rhea-id', 'annotation score', 'comment(CAUTION)', 'comment(MISCELLANEOUS)',
                  'keywords', 'context', 'existence', 'tools', 'reviewed', 'comment(BIOTECHNOLOGY)', 'go',
                  'go(biological process)', 'go(cellular component)', 'go(molecular function)', 'go-id',
                  'citationmapping', 'citation', 'families', 'database(Araport)', 'database(BioCyc)',
                  'database(BRENDA)', 'database(CDD)', 'database(CGD)', 'database(CollecTF)', 'database(dictyBase)',
                  'database(DNASU)', 'database(DrugBank)', 'database(DrugCentral)', 'database(EchoBASE)',
                  'database(eggNOG)', 'database(EMBL)', 'database(Ensembl)', 'database(EnsemblBacteria)',
                  'database(EnsemblFungi)', 'database(EnsemblPlants)', 'database(GeneID)', 'database(KEGG)',
                  'database(OrthoDB)', 'database(PATRIC)',
                  'database(Pfam)', 'database(PRO)',
                  'database(RefSeq)', 'database(SGD)',
                  'database(TIGRFAMs)', 'database(UniPathway)']


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


chunks_lof = []
submission_chunks = list(divide_chunks(random_saccs, submission_chunk_size))
for chunk in submission_chunks:
    chunk_frame = get_uniprot_sequences(uniprot_ids=submission_chunks[0], uniprot_cols=include_fields)
    chunks_lof.append(chunk_frame)

uniprot_frame = pd.concat(chunks_lof)

uniprot_frame.to_csv("target/uniprot_response_frame.tsv", sep="\t", index=False)

# # dtypedict or scalar, optional
# # Specifying the datatype for columns. If a dictionary is used, the keys should be the column names
# # and the values should be the SQLAlchemy types or strings for the sqlite3 legacy mode.
# # If a scalar is provided, it will be applied to all columns.
#
# # if_exists{‘fail’, ‘replace’, ‘append’}, default ‘fail’

uniprot_frame.to_sql(name="uniprot_annotations", con=seq2ids_db_conn, if_exists="replace", index=True)

# # ---
# # for record-keeping
#
# skip_fields = ['Virus hosts', 'database(CCDS)', 'database(BioGRID)', 'comment(PATHWAY)', 'comment(CATALYTIC ACTIVITY)',
#                'chebi', 'chebi(Catalytic activity)', 'chebi-id', 'features', 'comment(DISEASE)',
#                'comment(DISRUPTION PHENOTYPE)', 'comment(PHARMACEUTICAL)', 'Cross-reference (UCSC)',
#                'database(PseudoCAP)', 'database(PharmGKB)',
#                'database(PomBase)', 'database(Micado)', 'database(MoonDB)',
#                'database(MoonProt)', 'database(GeneReviews)', 'database(GO)', 'database(Gramene)',
#                'database(GuidetoPHARMACOLOGY)', 'database(ENZYME)', 'database(FlyBase)',
#                'database(GenBank)', "Erroneous termination", "Erroneous translation", "Frameshift", "Mass spectrometry",
#                "Polymorphism", "RNA editing"]
