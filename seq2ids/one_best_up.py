import pandas as pd
import sqlite3

pd.set_option("display.max_columns", None)

sqlite_file = "target/seq2ids.db"

con = sqlite3.connect(sqlite_file)
df = pd.read_sql_query(f"SELECT * from blast_results where blast_db = 'swissprot'", con)

qacc_max_bitscore = df.groupby("qacc")["bitscore"].max()

qacc_max_bitscore_frame = qacc_max_bitscore.to_frame()

qacc_max_bitscore_frame["qacc"] = qacc_max_bitscore_frame.index

qacc_max_bitscore_frame.reset_index(drop=True, inplace=True)

has_max_score = qacc_max_bitscore_frame.merge(right=df, how="inner")

qacc_counts = has_max_score["qacc"].value_counts()

one_best_qaccs = qacc_counts.loc[qacc_counts == 1]

one_best_frame = has_max_score.loc[
    has_max_score["qacc"].isin(list(one_best_qaccs.index))
]


for_insertion_frame = pd.DataFrame(one_best_frame["sacc"])
for_insertion_frame["minpos"] = None
for_insertion_frame["maxpos"] = None

for_insertion_frame.to_sql("one_best_up", con=con, if_exists="replace", index=False)