import datetime
import logging
import time

import click
import click_log
from Bio.Blast import NCBIWWW

from seq2ids.get_seqs_from_db import SeqsFromDb

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--secrets_file", type=click.Path(exists=True), required=True)
def seq2ids(secrets_file: str):
    """
    Gets slots, listed in config_tsv, from source_model and puts them in recipient_model
    :param secrets_file:
    :return:
    """

    # todo help for each option
    # todo better docstring

    sfd = SeqsFromDb()
    sfd.get_secrets_dict(secrets_file=secrets_file)
    rf = sfd.query_to_frame()
    # print(rf)

    seqs = rf['sequence'].tolist()

    first_seq = seqs[0]

    print(first_seq)

    # sequence_file = "blast_example.fasta"
    #
    # # make safer with with
    # sequence_data = open(sequence_file).read()
    #
    # print(sequence_data)

    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print(st)

    # seq_record = next(SeqIO.parse(sequence_file, 'fasta'))
    #
    # print(seq_record.id)
    # print(seq_record.seq)

    # one minute?
    result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=first_seq)

    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print(st)

    blast_results = result_handle.read()

    print(blast_results)


if __name__ == "__main__":
    seq2ids()
