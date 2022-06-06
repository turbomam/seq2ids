import logging
import pandas as pd

import click
import click_log

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click.option(
    "--input",
    default=None,
    type=click.Path(),
    required=True,
    help="Path to input TSV file",
)
@click.option(
    "--output",
    default=None,
    type=click.Path(),
    required=True,
    help="Path to output TSV file",
)
@click.option(
    "--annotation-cols",
    multiple=True,
    default=[
        "gene_name_1ry",
        "gene_names",
        "organism",
        "category",
        "prot_names",
        "prot_function",
        "all_go",
    ],
    help="list of annotation columns",
)
@click_log.simple_verbosity_option(logger)
def parts_best_match(input, output, annotation_cols):
    # read in spreadsheet which we need to trim to one blast hit
    # per seq_name
    df = pd.read_csv(input, header=0, sep="\t")

    # remove whole duplicates from original tsv
    df = df[df.nunique(1) > 1]

    # log whole duplicates if any
    # whole duplicates in this context are those rows
    # which are similar in all columns
    logger.debug(f"Whole duplicates in {input}:")
    logger.debug(df[df.nunique(1) <= 1])

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
    annotations_col_list = list(annotation_cols)

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
    df = df.drop(columns=["annotations_col_count", "weighted_score"])

    # write bitscore and presence of metadata filtered part_characterization.tsv
    # file into another tsv file in target folder
    df.to_csv(output, sep="\t", index=False)
    logger.info(f"Filtered {input} file written to: {output}")
