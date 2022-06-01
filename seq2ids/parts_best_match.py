import pandas as pd

# read in spreadsheet which we need to trim to one blast hit
# per seq_name
df = pd.read_csv("target/part_characterization.tsv", header=0, sep="\t")

# remove whole duplicates from original tsv
df = df[df.nunique(1) > 1]

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

# count number of filled columns for all bitscore filtered groups
df["count"] = df.notnull().sum(axis=1)

# pick seq_name row with highest count value
# i.e., the row which is richest in metadata
df = df.loc[df.reset_index().groupby(["seq_name"])["count"].idxmax()]

# drop the count column since it does not add any value to the
# data itself
df = df.drop(columns="count")

# write bitscore and presence of metadata filtered part_characterization.tsv
# file into another tsv file in target folder
df.to_csv("target/filtered_part_characterization.tsv", sep="\t")
