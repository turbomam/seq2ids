# todo help for each option
# todo better docstring
# todo combine flanks from any one given IF?

import logging

import click
import click_log
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seq2ids.get_seqs_from_db import SeqsFromDb

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
@click.option("--secrets_file", type=click.Path(exists=True), required=True)
@click.option("--fasta_out", type=click.Path(), required=True)
@click.option("--metadata_tsv_out", type=click.Path(), required=True)
@click.option("--min_len", type=int, help="inclusive")
@click.option("--max_len", type=int, help="exclusive")
def seq2ids(secrets_file: str, fasta_out: str, metadata_tsv_out: str, min_len: int, max_len: int):
    """
    Gets slots, listed in config_tsv, from source_model and puts them in recipient_model
    :param secrets_file:
    :param fasta_out:
    :param metadata_tsv_out:
    :param min_len:
    :param max_len:
    :return:
    """

    sfd = SeqsFromDb()
    sfd.get_secrets_dict(secrets_file=secrets_file)
    rf = sfd.query_to_frame()

    for_fasta = rf[['id', 'sequence']].copy()

    metadata = for_fasta.copy()
    metadata['seq_len'] = metadata['sequence'].str.len()
    metadata = metadata[['id', 'seq_len']]

    if min_len is None:
        min_len = 0
    if max_len is None:
        max_len = metadata['seq_len'].max()

    print(f"min {min_len}")
    print(f"max {max_len}")

    # was using seq_name
    desired_ids = metadata.loc[metadata['seq_len'].ge(min_len) & metadata['seq_len'].lt(max_len), 'id'].tolist()

    desired_seqs = for_fasta.loc[for_fasta['id'].isin(desired_ids)]

    # desired_seqs.sort_values('id')

    metadata.to_csv(metadata_tsv_out, sep='\t', index=False)

    ds_lod = desired_seqs.to_dict(orient='records')

    with open(fasta_out, 'w') as f_out:
        for seqs in ds_lod:
            sr = SeqRecord(Seq(seqs['sequence']), str(seqs['id']), '', '')
            r = SeqIO.write(sr, f_out, 'fasta')
            if r != 1:
                print('Error while writing sequence:  ' + sr.id)


if __name__ == "__main__":
    seq2ids()
