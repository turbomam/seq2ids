.PHONY: all clean

clean:
	rm -rf seq2ids.*

# 		--max_len 200 \
# sortin on ids into fasta not working
seq2ids.fasta:
	poetry run seq2ids \
		--secrets_file local/secrets.yaml \
		--fasta_out $@ \
		--min_len 30 \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)

seq2ids_under_30.fasta:
	poetry run seq2ids \
		--secrets_file local/secrets.yaml \
		--fasta_out $@ \
		--max_len 30 \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)


