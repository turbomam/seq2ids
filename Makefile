.PHONY: all clean blast_res_to_sqlite

clean:
	rm -rf seq2ids.*
	rm -rf target/*

# 		--max_len 200 \
# sorting on ids into fasta not working
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

# todo run blast

# 1 get sequences, names, part ids, types, lengths from felix postgres db
target/parts_sequences_plus.tsv:
	psql -h localhost -p 1111 -d felix -U mam -f sql/parts_sequences_plus.sql -F'	' --no-align --pset footer > $@
	sqlite3 target/seq2ids.db ".mode tabs" ".import target/parts_sequences_plus.tsv parts_sequences_plus" ""

# select
  #	"type" ,
  #	count(1)
  #from
  #	parts_sequences_plus psp
  #group by
  #	"type"
  #order by
  #	COUNT(1) desc ;

#|type     |count(1)|
#|---------|--------|
#|insertion|228     |
#|sequence |173     |
#|flank2   |97      |
#|deletion |97      |
#|flank1   |95      |


# todo long and short blast results!
# todo make sure they were run with the right blast tasks
# todo could programmatically generate .ini file's -outfmt and use the same for the initial header import below
# todo create VIEWs?
blast_res_to_sqlite: target/parts_sequences_plus.tsv
	sqlite3 target/seq2ids.db ".mode tabs" ".import seq2ids_elastic-blast_header.tsv blast_results" ""
	sqlite3 target/seq2ids.db ".mode tabs" ".import local/batch_000-blastn-nt.out_over.tsv  blast_results" ""
	sqlite3 target/seq2ids.db ".mode tabs" ".import local/batch_000-blastn-nt.out_under.tsv blast_results" ""
	sqlite3 target/seq2ids.db < sql/indices.sql
	# how many HSPs use each genome?
	sqlite3 target/seq2ids.db < sql/insertion_genome_hsp_counts.sql
	# rank each genome for each query
	sqlite3 target/seq2ids.db < sql/insertions_querys_genome_ranking.sql
	sqlite3 target/seq2ids.db < sql/insertion_genome_list.sql
	sqlite3 target/seq2ids.db < sql/ranges_to_download.sql


insertions_ids: clean target/parts_sequences_plus.tsv blast_res_to_sqlite target/insertion_genome_list.tsv

#wget -O /path/to/your.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=<acc[.ver]>"
