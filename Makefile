# doesn't yet do any special handling of deletion flanks
# may require ssh tunel to postgres

max_eval=1e-20
blast_thread_count=10

.PHONY: uniprot_approach load_seqs_blast_result clean blast_res_to_sqlite  nt_approach all

target/parts_partial.tsv:
	# may contain carriage returns
	psql -h localhost -p 1111 -d felix -U mam -f sql/parts_partial.sql -F'	' --no-align --pset footer > $@
	sqlite3 target/seq2ids.db ".mode tabs" ".import target/parts_partial.tsv parts_partial" ""

target/modifications.tsv:
	# may contain carriage returns
	psql -h localhost -p 1111 -d felix -U mam -f sql/modifications.sql -F'	' --no-align --pset footer > $@
	sqlite3 target/seq2ids.db ".mode tabs" ".import target/modifications.tsv modifications" ""

target/part_characterization.tsv: uniprot_approach
	sqlite3 target/seq2ids.db < sql/part_characterization.sql > $@

uniprot_approach: clean target/seq2ids_v_uniprot.tsv load_seqs_blast_result target/parts_partial.tsv target/modifications.tsv
	poetry run python seq2ids/get_uniprot_entries.py
	sqlite3 target/seq2ids.db < sql/parts_up_annotations.sql

clean:
	#rm -rf seq2ids.*
	rm -rf target/*
	rm -rf blastdbs/*
	rm -rf swissprot*
	rm -rf taxdb*
	#rm -rf seq2ids_v_uniprot.tsv

# run blast
#   locally?
#   submissions against nt, to NCBI, via BioPython SLOW
#   elastic blast $, not instantaneous either

  # swissprot-prot-metadata.json               2022-04-04 15:14  444
# swissprot.tar.gz                           2022-04-04 15:14  180M
  # swissprot.tar.gz.md5                       2022-04-04 15:14   51
  # taxdb-metadata.json                        2021-06-07 10:40  388
# taxdb.tar.gz

blastdbs/swissprot.psq:
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
	tar -xzvf swissprot.tar.gz --directory blastdbs

taxdb.bti: blastdbs/swissprot.psq
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
	tar -xzvf taxdb.tar.gz --directory .

# assumes blastx is on the path
# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
target/seq2ids_v_uniprot.tsv: target/seq2ids.fasta taxdb.bti
	blastx \
		-query $< \
		-db blastdbs/swissprot \
		-num_threads ${blast_thread_count} \
		-out target/at_delim_blast.txt \
		-evalue ${max_eval} \
		-outfmt "6 delim=@ qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles staxids sskingdoms sscinames sblastnames scomnames"
	cat target/at_delim_blast.txt | tr '@' '\t' > $@

# evalue bitscore length pident qacc qcovhsp qcovs qstart qend qseqid sacc sallacc sallgi sallseqid salltitles sblastnames scomnames sstart send sgi sseqid staxids stitle sskingdoms ssciname gapopen mismatch"
# qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles staxids sskingdoms sscinames sblastnames scomnames"
# diff: gapopen mismatch qseqid sallgi score sgi

# sorting on ids into fasta not working
target/seq2ids.fasta:
	poetry run seq2ids \
		--secrets_file local/secrets.yaml \
		--fasta_out $@ \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)

load_seqs_blast_result: target/parts_sequences_plus.tsv blast_res_to_sqlite

# 1 get sequences, names, part ids, types, lengths from felix postgres db
target/parts_sequences_plus.tsv:
	psql -h localhost -p 1111 -d felix -U mam -f sql/parts_sequences_plus.sql -F'	' --no-align --pset footer > $@
	sqlite3 target/seq2ids.db ".mode tabs" ".import target/parts_sequences_plus.tsv parts_sequences_plus" ""

# todo long and short blast results!
# todo make sure they were run with the right blast tasks
# todo could programmatically generate .ini file's -outfmt and use the same for the initial header import below
# todo create VIEWs?
# todo when ranking, is the precedence ltr?
blast_res_to_sqlite: target/parts_sequences_plus.tsv
	#sqlite3 target/seq2ids.db ".mode tabs" ".import seq2ids_elastic-blast_header.tsv blast_results" ""
	sqlite3 target/seq2ids.db < sql/create_blast_results_table.sql
##	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_repeat/batch_000-blastn-nt.out  blast_results" ""
##	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_under_30_repeat/batch_000-blastn-nt.out blast_results" ""
##	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_all_swissprot/batch_000-blastx-swissprot.out blast_results" ""
##	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_all_swissprot/batch_001-blastx-swissprot.out blast_results" ""
##	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_all_swissprot/batch_002-blastx-swissprot.out blast_results" ""
#	sqlite3 target/seq2ids.db ".mode tabs" ".import local/seq2ids_all_swissprot/batch_002-blastx-swissprot.out blast_results" ""
	sqlite3 target/seq2ids.db ".mode tabs" ".import target/seq2ids_v_uniprot.tsv blast_results" ""
	sqlite3 target/seq2ids.db < sql/indices.sql
	# how many HSPs use each genome?
	sqlite3 target/seq2ids.db < sql/insertion_genome_hsp_counts.sql
	# rank each genome for each query
	sqlite3 target/seq2ids.db < sql/insertions_querys_genome_ranking.sql
	sqlite3 target/seq2ids.db < sql/ranges_to_download.sql

target/parts_partial.tsv:
	# may contain carriage returns
	psql -h localhost -p 1111 -d felix -U mam -f sql/parts_partial.sql -F'	' --no-align --pset footer > $@
	sqlite3 target/seq2ids.db ".mode tabs" ".import target/parts_partial.tsv parts_partial" ""

target/modifications.tsv:
	# may contain carriage returns
	psql -h localhost -p 1111 -d felix -U mam -f sql/modifications.sql -F'	' --no-align --pset footer > $@
	sqlite3 target/seq2ids.db ".mode tabs" ".import target/modifications.tsv modifications" ""

# ---

seq2ids_under_30_nt.fasta:
	poetry run seq2ids \
		--secrets_file local/secrets.yaml \
		--fasta_out $@ \
		--max_len 30 \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)

seq2ids_30_nt_plus.fasta:
	poetry run seq2ids \
		--secrets_file local/secrets.yaml \
		--fasta_out $@ \
		--max_len 30 \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)

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

nt_approach: load_seqs_blast_result
	sqlite3 target/seq2ids.db < sql/smin_smax.sql
	poetry run python seq2ids/efetch_features.py
	sqlite3 target/seq2ids.db < sql/b2f_summary.sql

# blastx (on path) installed from...

