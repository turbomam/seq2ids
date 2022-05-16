import logging

import click
import click_log

import pandas as pd

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--tsv_in", type=click.Path(exists=True), required=True)
@click.option("--tsv_out", type=click.Path(), required=True)
@click.option("--col_val", required=True)
def cli(tsv_in: str, tsv_out: str, col_val: str):
    """
    Add a column with a constant value to a TSV file
    :param tsv_in:
    :param tsv_out:
    :param col_name:
    :param col_val:
    :return:
    """

    in_frame = pd.read_csv(tsv_in, sep="\t", header=None)
    in_frame["added"] = col_val
    in_frame.to_csv(tsv_out, sep="\t", index=False, header=False)

    pass


if __name__ == "__main__":
    cli()
