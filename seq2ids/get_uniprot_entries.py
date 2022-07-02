# import pprint
import sqlite3

import pandas as pd
import requests

pd.set_option("display.max_columns", None)

seq2ids_db_fp = "target/seq2ids.db"

paname = "primaryAccession"

# todo better handling of fpbase esp annotations and prioritizing relative to uniprot
# todo better handling of coimpound insertions
# todo ready to start including other non-flank sequence types?


def get_highlights(obj):
    highlights = {}

    highlights["rec_full_name"] = obj["proteinDescription"]["recommendedName"][
        "fullName"
    ]["value"]
    # todo more names symbols synonyms etc.?
    if "genes" in obj:
        for current_gene_info in obj["genes"]:
            if "geneName" in current_gene_info:
                # print(yaml.dump(current_gene_info['geneName']['value']))
                highlights["geneName"] = current_gene_info["geneName"]["value"]

    if "ecNumbers" in obj["proteinDescription"]["recommendedName"]:
        ec_temp = [
            i["value"]
            for i in obj["proteinDescription"]["recommendedName"]["ecNumbers"]
        ]
        ec_list = "|".join(ec_temp)
        highlights["ecNumbers"] = ec_list

    go_ids = []
    go_labels = []
    for xref in obj["uniProtKBCrossReferences"]:
        if xref["database"] == "GO":
            go_ids.append(xref["id"])
            # print(xref['id'])
            for prop in xref["properties"]:
                if prop["key"] == "GoTerm":
                    go_labels.append(prop["value"])

    if go_ids:
        highlights["go_ids"] = "|".join(go_ids)
        highlights["go_labels"] = "|".join(go_labels)

    highlights["uniProtkbId"] = obj["uniProtkbId"]
    highlights["annotationScore"] = obj["annotationScore"]
    highlights["primaryAccession"] = obj["primaryAccession"]
    highlights["org_scientificName"] = obj["organism"]["scientificName"]
    highlights["taxonId"] = obj["organism"]["taxonId"]

    return highlights


seq2ids_db_conn = sqlite3.connect(seq2ids_db_fp)

annotatable_q = f"select distinct sacc from blast_results where pident > 90 and qcovs > 90 and blast_db = 'swissprot'"

annotatable_res = pd.read_sql_query(annotatable_q, seq2ids_db_conn)

saccs = annotatable_res["sacc"].to_list()
saccs.sort()

annotations = {}
for sacc in saccs:
    print(sacc)
    url = f"https://rest.uniprot.org/uniprotkb/{sacc}.json"
    sacc_res = requests.get(url).json()
    annotations[sacc] = sacc_res

annotation_ids = list(annotations.keys())
annotation_ids.sort()

highlight_list = []
for current_id in annotation_ids:
    print(current_id)
    highlight = get_highlights(annotations[current_id])
    highlight_list.append(highlight)

highlight_frame = pd.DataFrame(highlight_list)

temp = highlight_frame.columns
temp = list(temp)

temp.remove(paname)
temp.sort()
temp = [paname] + temp

highlight_frame = highlight_frame[temp]

seq2ids_db_conn = sqlite3.connect(seq2ids_db_fp)

highlight_frame.to_sql(
    name="uniprot_annotations", con=seq2ids_db_conn, if_exists="replace", index=False
)


# include_fields = ['id', 'entry name', 'genes', 'genes(PREFERRED)', 'genes(ALTERNATIVE)', 'genes(OLN)', 'genes(ORF)',
#                   'organism-id', 'organism', 'lineage(ALL)', 'protein names', 'proteome', 'database(BioCyc)',
#                   'database(BRENDA)', 'database(CDD)', 'comment(ALTERNATIVE PRODUCTS)',
#                   'comment(ERRONEOUS GENE MODEL PREDICTION)', 'comment(ERRONEOUS INITIATION)',
#                   'comment(SEQUENCE CAUTION)', 'comment(ENZYME REGULATION)', 'comment(FUNCTION)',
#                   'ec', 'rhea-id', 'annotation score', 'comment(CAUTION)', 'comment(MISCELLANEOUS)',
#                   'keywords', 'context', 'existence', 'tools', 'reviewed', 'comment(BIOTECHNOLOGY)', 'go',
#                   'go(biological process)', 'go(cellular component)', 'go(molecular function)', 'go-id',
#                   'citationmapping', 'citation', 'families', 'database(Araport)', 'database(BioCyc)',
#                   'database(BRENDA)', 'database(CDD)', 'database(CGD)', 'database(CollecTF)', 'database(dictyBase)',
#                   'database(DNASU)', 'database(DrugBank)', 'database(DrugCentral)', 'database(EchoBASE)',
#                   'database(eggNOG)', 'database(EMBL)', 'database(Ensembl)', 'database(EnsemblBacteria)',
#                   'database(EnsemblFungi)', 'database(EnsemblPlants)', 'database(GeneID)', 'database(KEGG)',
#                   'database(OrthoDB)', 'database(PATRIC)',
#                   'database(Pfam)', 'database(PRO)',
#                   'database(RefSeq)', 'database(SGD)',
#                   'database(TIGRFAMs)', 'database(UniPathway)']


# def get_uniprot_sequences(uniprot_ids: List, uniprot_cols: List) -> pd.DataFrame:
#     """
#     https://www.biostars.org/p/94422/
#
#     Retrieve uniprot sequences based on a list of uniprot sequence identifier.
#
#     For large lists it is recommended to perform batch retrieval.
#
#     documentation which columns are available:
#     https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames
#
#     this python script is based on
#     https://www.biostars.org/p/67822/
#
#     Parameters:
#         uniprot_ids: List, list of uniprot identifier
#         uniprot_cols: List, list of display columns (https://www.uniprot.org/help/uniprotkb_column_names and https://www.uniprot.org/docs/dbxref)
#
#     Returns:
#         pd.DataFrame, pandas dataframe with uniprot id column and sequence
#     """
#
#     print(uniprot_ids)
#
#     print(uniprot_cols)
#
#     joined_cols = ",".join(uniprot_cols)
#
#     url = "https://www.uniprot.org/uploadlists/"  # This is the webserver to retrieve the Uniprot data
#     params = {
#         "from": "ACC",
#         "to": "ACC",
#         "format": "tab",
#         "query": " ".join(uniprot_ids),
#         "columns": joined_cols,
#     }
#
#     response = requests.get(url, params=params)
#
#     print(response.content)
#
#     # b'<html>\r\n<head><title>404 Not Found</title></head>\r\n<body>\r\n<center><h1>404 Not Found</h1></center>\r\n<hr><center>nginx/1.21.6</center>\r\n</body>\r\n</html>\r\n'
#
#     # print(response.request.url)
#
#     response_frame = pd.read_csv(io.BytesIO(response.content), sep="\t")
#     response_frame = response_frame.iloc[:, :-1]
#
#     return response_frame
#
#
# def divide_chunks(l, n):
#     # looping till length l
#     for i in range(0, len(l), n):
#         yield l[i : i + n]


# # # # # dtypedict or scalar, optional
# # # # # Specifying the datatype for columns. If a dictionary is used, the keys should be the column names
# # # # # and the values should be the SQLAlchemy types or strings for the sqlite3 legacy mode.
# # # # # If a scalar is provided, it will be applied to all columns.
# # # #
# # # # # if_exists{‘fail’, ‘replace’, ‘append’}, default ‘fail’
# # # # # for record-keeping
# # # #
# # # # skip_fields = ['Virus hosts', 'database(CCDS)', 'database(BioGRID)', 'comment(PATHWAY)', 'comment(CATALYTIC ACTIVITY)',
# # # #                'chebi', 'chebi(Catalytic activity)', 'chebi-id', 'features', 'comment(DISEASE)',
# # # #                'comment(DISRUPTION PHENOTYPE)', 'comment(PHARMACEUTICAL)', 'Cross-reference (UCSC)',
# # # #                'database(PseudoCAP)', 'database(PharmGKB)',
# # # #                'database(PomBase)', 'database(Micado)', 'database(MoonDB)',
# # # #                'database(MoonProt)', 'database(GeneReviews)', 'database(GO)', 'database(Gramene)',
# # # #                'database(GuidetoPHARMACOLOGY)', 'database(ENZYME)', 'database(FlyBase)',
# # # #                'database(GenBank)', "Erroneous termination", "Erroneous translation", "Frameshift", "Mass spectrometry",
# # # #                "Polymorphism", "RNA editing"]
