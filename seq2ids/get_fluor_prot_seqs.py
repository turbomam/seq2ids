import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

fp_url = "https://www.fpbase.org/api/proteins/?format=json"
fasta_out = "data/fpbase.fasta"

fp = requests.get(fp_url)

fp_lod = fp.json()

with open(fasta_out, "w") as f_out:
    for seqs in fp_lod:
        if seqs['slug'] and seqs['seq']:
            tidy = re.sub(r'\s+', '', seqs["seq"])
            sr = SeqRecord(Seq(tidy), str(seqs["slug"]), "", "")
            r = SeqIO.write(sr, f_out, "fasta")
            if r != 1:
                print("Error while writing sequence:  " + sr.id)
