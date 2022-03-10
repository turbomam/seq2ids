from Bio.Blast import NCBIWWW
from Bio import SeqIO

import logging

from typing import Optional, Dict, List, Any

import click
import click_log

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--word")
# @click.option("--right_model", type=click.Path(exists=True), required=True)
# @click.option("--yaml_output", type=click.Path(), required=True)
def seq2ids(word: str):
    """
    Gets slots, listed in config_tsv, from source_model and puts them in recipient_model
    :param word:
    :return:
    """

    # todo help for each option
    # todo docstring

    print(word)

    sequence_file = "blast_example.fasta"

    # make safer with with
    sequence_data = open(sequence_file).read()

    print(sequence_data)

    # seq_record = next(SeqIO.parse(sequence_file, 'fasta'))
    #
    # print(seq_record.id)
    # print(seq_record.seq)

    result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)


if __name__ == "__main__":
    seq2ids()
