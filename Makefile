.PHONY: all clean blast_res_to_sqlite load_seqs_blast_result nt_approach uniprot_approach

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
# todo when renking, is the precedence ltr?
blast_res_to_sqlite: target/parts_sequences_plus.tsv
	#sqlite3 target/seq2ids.db ".mode tabs" ".import seq2ids_elastic-blast_header.tsv blast_results" ""
	sqlite3 target/seq2ids.db < sql/create_blast_results_table.sql
#	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_repeat/batch_000-blastn-nt.out  blast_results" ""
#	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_under_30_repeat/batch_000-blastn-nt.out blast_results" ""
	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_all_swissprot/batch_000-blastx-swissprot.out blast_results" ""
	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_all_swissprot/batch_001-blastx-swissprot.out blast_results" ""
	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_all_swissprot/batch_002-blastx-swissprot.out blast_results" ""
	sqlite3 target/seq2ids.db < sql/indices.sql
	# how many HSPs use each genome?
	sqlite3 target/seq2ids.db < sql/insertion_genome_hsp_counts.sql
	# rank each genome for each query
	sqlite3 target/seq2ids.db < sql/insertions_querys_genome_ranking.sql
	sqlite3 target/seq2ids.db < sql/ranges_to_download.sql


load_seqs_blast_result: clean target/parts_sequences_plus.tsv blast_res_to_sqlite

nt_approach: load_seqs_blast_result
	sqlite3 target/seq2ids.db < sql/smin_smax.sql
	poetry run python seq2ids/efetch_features.py
	sqlite3 target/seq2ids.db < sql/b2f_summary.sql


uniprot_approach: load_seqs_blast_result
	poetry run python seq2ids/get_uniprot_entries.py
	sqlite3 target/seq2ids.db < sql/parts_up_annotations.sql




