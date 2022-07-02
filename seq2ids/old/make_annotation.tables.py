# import pprint
import sqlite3
from pathlib import Path

import pandas as pd
import yaml

seq2ids_db_fp = "target/seq2ids.db"

paname = "primaryAccession"

annotations_yaml = "up.yaml"

annotations = yaml.safe_load(Path(annotations_yaml).read_text())

annotation_ids = list(annotations.keys())
annotation_ids.sort()


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
