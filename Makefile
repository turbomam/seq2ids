# doesn't yet do any special handling of deletion flanks
# may require ssh tunel to postgres

max_eval=1e-20
blast_thread_count=10
#selected_sqlite_input_db=local/felix_dump.db # prod for your eyes only
selected_sqlite_input_db=data/sqlite_not_postgres.db # test for public
destination_sqlite_db=target/seq2ids.db

.PHONY: uniprot_approach load_seqs_blast_result clean blast_res_to_sqlite  nt_approach all sqlite_input live_db

all: target/part_characterization.tsv

target/part_characterization.tsv: uniprot_sqlite_input
	sqlite3 $(destination_sqlite_db) < sql/part_characterization.sql > $@

uniprot_sqlite_input: clean local/felix_dump.db live_db target/seq2ids_v_uniprot.tsv target/seq2ids_v_fpbase.tsv blast_res_to_sqlite
	poetry run python seq2ids/get_uniprot_entries.py
	sqlite3 $(destination_sqlite_db) < sql/parts_up_annotations.sql

live_db:
	cp $(selected_sqlite_input_db) local/live_sqlite.db

clean:
	rm -rf target/*
	rm -rf local/felix_dump.db

squeaky_clean: clean
	#rm -rf seq2ids.*
	rm -rf blastdbs/*
	rm -rf swissprot*
	rm -rf taxdb*
	#	rm -rf seq2ids_v_uniprot.tsv

blastdbs/swissprot.psq:
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
	tar -xzvf swissprot.tar.gz --directory blastdbs

taxdb.bti: blastdbs/swissprot.psq
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
	tar -xzvf taxdb.tar.gz --directory .

 # FOR SQLITE INPUT
target/seq2ids.fasta:
	poetry run seq2ids \
		--sqlite_file local/live_sqlite.db \
		--fasta_out $@ \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)

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
	cat target/at_delim_blast.txt | tr '@' '\t' > target/tab_delim_blast.tsv
	poetry run python seq2ids/add_col.py \
		--tsv_in target/tab_delim_blast.tsv \
		--tsv_out $@ \
		--col_val swissprot
	rm target/at_delim_blast.txt target/tab_delim_blast.tsv

target/seq2ids_v_fpbase.tsv: target/seq2ids.fasta data/fpbase.fasta.psq
	blastx \
		-query $< \
		-db data/fpbase.fasta \
		-num_threads ${blast_thread_count} \
		-out target/at_delim_blast.txt \
		-evalue ${max_eval} \
		-outfmt "6 delim=@ qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles staxids sskingdoms sscinames sblastnames scomnames"
	cat target/at_delim_blast.txt | tr '@' '\t' > target/tab_delim_blast.tsv
	poetry run python seq2ids/add_col.py \
		--tsv_in target/tab_delim_blast.tsv \
		--tsv_out $@ \
		--col_val fpbase
	rm target/at_delim_blast.txt target/tab_delim_blast.tsv

# todo separate nt blast searches for > 30 nt vs shorter queries
# todo make sure they were run with the right blast subtasks
# todo could programmatically generate elastic blast's .ini file's -outfmt
#   and use the same for the initial header import below
# todo create VIEWs?
# todo when ranking, is the precedence ltr?
blast_res_to_sqlite: target/seq2ids_v_uniprot.tsv target/seq2ids_v_fpbase.tsv
	#sqlite3 $(destination_sqlite_db) ".mode tabs" ".import seq2ids_elastic-blast_header.tsv blast_results" ""
	sqlite3 $(destination_sqlite_db) < sql/create_blast_results_table.sql
## this is an example of parsing elastic blast results
#	sqlite3 $(destination_sqlite_db) ".mode tabs" ".import local/seq2ids_repeat/batch_000-blastn-nt.out  blast_results" ""
	sqlite3 $(destination_sqlite_db) ".mode tabs" ".import target/seq2ids_v_uniprot.tsv blast_results" ""
	sqlite3 $(destination_sqlite_db) ".mode tabs" ".import target/seq2ids_v_fpbase.tsv blast_results" ""
	sqlite3 $(destination_sqlite_db) < sql/by_attachement.sql
	sqlite3 $(destination_sqlite_db) < sql/indices.sql
	# how many HSPs use each genome?
	sqlite3 $(destination_sqlite_db) < sql/insertion_genome_hsp_counts.sql
	# rank each genome for each query
	sqlite3 $(destination_sqlite_db) < sql/insertions_querys_genome_ranking.sql
	sqlite3 $(destination_sqlite_db) < sql/ranges_to_download.sql
	# seq2ids/one_best_up.py undoes several previous steps that try to identify the minimal number of uniprot annotations that need to be retrieved
	poetry run python seq2ids/one_best_up.py

data/fpbase.fasta:
	poetry run python seq2ids/get_fluor_prot_seqs.py

data/fpbase.fasta.psq: data/fpbase.fasta
	makeblastdb -in $< -dbtype prot

local/felix_dump.db:
	poetry run sh bash/pgsql2sqlite.sh mam 1111 parts,parts_sequences,modifications

# ----


target/parts_partial.tsv:
	# may contain carriage returns
	psql -h localhost -p 1111 -d felix -U mam -f sql/parts_partial.sql -F'	' --no-align --pset footer > $@
	sqlite3 $(destination_sqlite_db) ".mode tabs" ".import target/parts_partial.tsv parts_partial" ""

target/modifications.tsv:
	# may contain carriage returns
	psql -h localhost -p 1111 -d felix -U mam -f sql/modifications.sql -F'	' --no-align --pset footer > $@
	sqlite3 $(destination_sqlite_db) ".mode tabs" ".import target/modifications.tsv modifications" ""

#target/part_characterization.tsv: uniprot_approach
#	sqlite3 $(destination_sqlite_db) < sql/part_characterization.sql > $@

uniprot_approach: clean target/seq2ids_v_uniprot.tsv load_seqs_blast_result target/parts_partial.tsv target/modifications.tsv
	poetry run python seq2ids/get_uniprot_entries.py
	sqlite3 $(destination_sqlite_db) < sql/parts_up_annotations.sql

# run blast
#   locally?
#   submissions against nt, to NCBI, via BioPython SLOW
#   elastic blast $, not instantaneous either

  # swissprot-prot-metadata.json               2022-04-04 15:14  444
# swissprot.tar.gz                           2022-04-04 15:14  180M
  # swissprot.tar.gz.md5                       2022-04-04 15:14   51
  # taxdb-metadata.json                        2021-06-07 10:40  388
# taxdb.tar.gz

# evalue bitscore length pident qacc qcovhsp qcovs qstart qend qseqid sacc sallacc sallgi sallseqid salltitles sblastnames scomnames sstart send sgi sseqid staxids stitle sskingdoms ssciname gapopen mismatch"
# qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles staxids sskingdoms sscinames sblastnames scomnames"
# diff: gapopen mismatch qseqid sallgi score sgi

## sorting on ids into fasta not working
## FOR POSTGRES INPUT
#target/seq2ids.fasta:
#	poetry run seq2ids \
#		--postgres_secrets_file local/secrets.yaml \
#		--fasta_out $@ \
#		--metadata_tsv_out $(subst .fasta,.tsv,$@)

# target/parts_sequences_plus.tsv required for Postgres approach
load_seqs_blast_result:  blast_res_to_sqlite

# 1 get sequences, names, part ids, types, lengths from felix postgres db
target/parts_sequences_plus.tsv:
	# # required for postgres apporach
	# psql -h localhost -p 1111 -d felix -U mam -f sql/parts_sequences_plus.sql -F'	' --no-align --pset footer > $@
	sqlite3 $(destination_sqlite_db) ".mode tabs" ".import target/parts_sequences_plus.tsv parts_sequences_plus" ""



# ---

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
	sqlite3 $(destination_sqlite_db) < sql/smin_smax.sql
	poetry run python seq2ids/efetch_features.py
	sqlite3 $(destination_sqlite_db) < sql/b2f_summary.sql

# blastx (on path) installed from...


# requires commenting out the target/seq2ids_v_uniprot.tsv version you DON'T want
postgres_input: clean target/seq2ids_v_uniprot.tsv

