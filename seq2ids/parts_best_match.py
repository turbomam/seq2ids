import pandas as pd

# read in spreadsheet which we need to trim to one blast hit
# per seq_name
df = pd.read_csv("target/part_characterization.tsv", header=0, sep="\t")

# remove whole duplicates from original tsv
df = df[df.nunique(1) > 1]

# drop rows where the bitscore value <= 0
df = df[df["bitscore"] > 0]

# sort the rows for each seq_name based on bitscore column
df = df.groupby("seq_name", as_index=False).apply(
    lambda x: x.sort_values("bitscore", ascending=False)
)

# select only the first two rows from each of the sorted groups
df = df.groupby("seq_name", as_index=False).apply(lambda x: x.head(2))

# applying groupby() multiple times results in a multi index
# being created which needs to be removed
df = df.reset_index()
df = df[df.columns.drop(list(df.filter(regex="level_")))]

# harcoded list of annotations columns
# TODO: this list should be read from a file input produced
# by the Uniprot annotations fetching script
annotations_col_list = [
    "gene_name_1ry",
    "gene_names",
    "organism",
    "category",
    "prot_names",
    "prot_function",
    "all_go",
]

# count number of filled columns for all bitscore filtered groups
df["annotations_col_count"] = df[annotations_col_list].notnull().sum(axis=1)

# compute weighted bitscore column count score
df["weighted_score"] = df["bitscore"] * (
    df["annotations_col_count"] / len(annotations_col_list)
)

# drop all rows where weighted_score = NaN
df.dropna(subset=["weighted_score"], how="all", inplace=True)

# pick seq_name row with highest weighted_score value
# i.e., the row which is richest in metadata
df = df.loc[df.groupby(["seq_name"])["weighted_score"].idxmax()]

# drop the annotations_col_count column since it does not add any value to the
# data itself
df = df.drop(columns="annotations_col_count")

# write bitscore and presence of metadata filtered part_characterization.tsv
# file into another tsv file in target folder
df.to_csv("target/filtered_part_characterization.tsv", sep="\t", index=False)
