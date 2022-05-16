# todo help for each option
# todo better docstring
# todo combine flanks from any one given IF?

import logging

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import click_log
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

# from get_seqs_from_db import SeqsFromDb

from seq2ids.get_seqs_from_db import SeqsFromDb

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command()
@click_log.simple_verbosity_option(logger)
# , required=False
# try https://pypi.org/project/click-option-group/
@optgroup.group('Input data sources', cls=RequiredMutuallyExclusiveOptionGroup,
                help='The sources of the input data')
@optgroup.option("--postgres_secrets_file", type=click.Path(exists=True))
@optgroup.option("--sqlite_file", type=click.Path(exists=True))
@click.option("--fasta_out", type=click.Path(), required=True)
@click.option("--metadata_tsv_out", type=click.Path(), required=True)
@click.option("--min_len", type=int, help="inclusive")
@click.option("--max_len", type=int, help="exclusive")
def seq2ids(postgres_secrets_file: str, sqlite_file: str, fasta_out: str, metadata_tsv_out: str, min_len: int,
            max_len: int):
    """
    Gets slots, listed in config_tsv, from source_model and puts them in recipient_model
    :param postgres_secrets_file:
    :param fasta_out:
    :param metadata_tsv_out:
    :param min_len:
    :param max_len:
    :return:
    """

    sfd = SeqsFromDb()

    if postgres_secrets_file:
        sfd.get_secrets_dict(secrets_file=postgres_secrets_file)
        rf = sfd.postgres_to_frame()
    elif sqlite_file:
        rf = sfd.sqlite_to_frame(sqlite_file)
    else:
        exit()

    for_fasta = rf[['id', 'sequence']].copy()
    for_fasta['sequence'] = for_fasta['sequence'].str.replace(' ', '')

    metadata = for_fasta.copy()
    metadata['seq_len'] = metadata['sequence'].str.len()
    metadata = metadata[['id', 'seq_len']]

    logger.info(f"sequence count        : {len(metadata.index)}")
    logger.info(f"min specified seq len : {min_len}")
    logger.info(f"min observed seq len  : {metadata['seq_len'].min()}")
    logger.info(f"max specified seq len : {max_len}")
    logger.info(f"max observed seq len  : {metadata['seq_len'].max()}")

    if not min_len:
        min_len = 0
    if not max_len:
        max_len = metadata['seq_len'].max() + 1

    # was using seq_name
    desired_ids = metadata.loc[metadata['seq_len'].ge(min_len) & metadata['seq_len'].lt(max_len), 'id'].tolist()

    desired_seqs = for_fasta.loc[for_fasta['id'].isin(desired_ids)]

    # desired_seqs.sort_values('id')

    metadata.to_csv(metadata_tsv_out, sep='\t', index=False)

    ds_lod = desired_seqs.to_dict(orient='records')

    with open(fasta_out, 'w') as f_out:
        for seqs in ds_lod:
            tidy = re.sub(r'\s+', '', seqs["sequence"])
            sr = SeqRecord(Seq(tidy), str(seqs['id']), '', '')
            r = SeqIO.write(sr, f_out, 'fasta')
            if r != 1:
                logger.error('Error while writing sequence:  ' + sr.id)


if __name__ == "__main__":
    seq2ids()
