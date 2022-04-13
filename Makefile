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
# todo when renking, is the precedence ltr?
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
# database from...

# blastx -query seq2ids.fasta -db /Users/MAM/swissprot/swissprot -num_threads 14 -out seq2ids_v_uniprot.tsv -evalue 1e-20 -outfmt "6 evalue bitscore length pident qacc qcovhsp qcovs qstart qend qseqid sacc sallacc sallgi sallseqid salltitles sblastnames scomnames sstart send sgi sseqid staxids stitle sskingdoms ssciname gapopen mismatch"

# -task options: blastx and

# blastx --help
  #USAGE
  #  blastx [-h] [-help] [-import_search_strategy filename]
  #    [-export_search_strategy filename] [-task task_name] [-db database_name]
  #    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
  #    [-negative_gilist filename] [-negative_seqidlist filename]
  #    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
  #    [-negative_taxidlist filename] [-ipglist filename]
  #    [-negative_ipglist filename] [-entrez_query entrez_query]
  #    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
  #    [-subject subject_input_file] [-subject_loc range] [-query input_file]
  #    [-out output_file] [-evalue evalue] [-word_size int_value]
  #    [-gapopen open_penalty] [-gapextend extend_penalty]
  #    [-qcov_hsp_perc float_value] [-max_hsps int_value]
  #    [-xdrop_ungap float_value] [-xdrop_gap float_value]
  #    [-xdrop_gap_final float_value] [-searchsp int_value]
  #    [-sum_stats bool_value] [-max_intron_length length] [-seg SEG_options]
  #    [-soft_masking soft_masking] [-matrix matrix_name]
  #    [-threshold float_value] [-culling_limit int_value]
  #    [-best_hit_overhang float_value] [-best_hit_score_edge float_value]
  #    [-subject_besthit] [-window_size int_value] [-ungapped] [-lcase_masking]
  #    [-query_loc range] [-strand strand] [-parse_deflines]
  #    [-query_gencode int_value] [-outfmt format] [-show_gis]
  #    [-num_descriptions int_value] [-num_alignments int_value]
  #    [-line_length line_length] [-html] [-sorthits sort_hits]
  #    [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
  #    [-num_threads int_value] [-mt_mode int_value] [-remote]
  #    [-comp_based_stats compo] [-use_sw_tback] [-version]

# DESCRIPTION
  #   Translated Query-Protein Subject BLAST 2.12.0+
  #
  #OPTIONAL ARGUMENTS
  # -h
  #   Print USAGE and DESCRIPTION;  ignore all other parameters
  # -help
  #   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
  # -version
  #   Print version number;  ignore other arguments
  #
  # *** Input query options
  # -query <File_In>
  #   Input file name
  #   Default = `-'
  # -query_loc <String>
  #   Location on the query sequence in 1-based offsets (Format: start-stop)
  # -strand <String, `both', `minus', `plus'>
  #   Query strand(s) to search against database/subject
  #   Default = `both'
  # -query_gencode <Integer, values between: 1-6, 9-16, 21-31, 33>
  #   Genetic code to use to translate query (see
  #   https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=
  #   cgencodes for details)
  #   Default = `1'
  #
  # *** General search options
  # -task <String, Permissible values: 'blastx' 'blastx-fast' >
  #   Task to execute
  #   Default = `blastx'
  # -db <String>
  #   BLAST database name
  #    * Incompatible with:  subject, subject_loc
  # -out <File_Out>
  #   Output file name
  #   Default = `-'
  # -evalue <Real>
  #   Expectation value (E) threshold for saving hits
  #   Default = `10'
  # -word_size <Integer, >=2>
  #   Word size for wordfinder algorithm
  # -gapopen <Integer>
  #   Cost to open a gap
  # -gapextend <Integer>
  #   Cost to extend a gap
  # -max_intron_length <Integer, >=0>
  #   Length of the largest intron allowed in a translated nucleotide sequence
  #   when linking multiple distinct alignments
  #   Default = `0'
  # -matrix <String>
  #   Scoring matrix name (normally BLOSUM62)
  # -threshold <Real, >=0>
  #   Minimum word score such that the word is added to the BLAST lookup table
  # -comp_based_stats <String>
  #   Use composition-based statistics:
  #       D or d: default (equivalent to 2 )
  #       0 or F or f: No composition-based statistics
  #       1: Composition-based statistics as in NAR 29:2994-3005, 2001
  #       2 or T or t : Composition-based score adjustment as in Bioinformatics
  #   21:902-911,
  #       2005, conditioned on sequence properties
  #       3: Composition-based score adjustment as in Bioinformatics 21:902-911,
  #       2005, unconditionally
  #   Default = `2'
  #
  # *** BLAST-2-Sequences options
  # -subject <File_In>
  #   Subject sequence(s) to search
  #    * Incompatible with:  db, gilist, seqidlist, negative_gilist,
  #   negative_seqidlist, taxids, taxidlist, negative_taxids, negative_taxidlist,
  #   ipglist, negative_ipglist, db_soft_mask, db_hard_mask
  # -subject_loc <String>
  #   Location on the subject sequence in 1-based offsets (Format: start-stop)
  #    * Incompatible with:  db, gilist, seqidlist, negative_gilist,
  #   negative_seqidlist, taxids, taxidlist, negative_taxids, negative_taxidlist,
  #   ipglist, negative_ipglist, db_soft_mask, db_hard_mask, remote
  #
  # *** Formatting options
  # -outfmt <String>
  #   alignment view options:
  #     0 = Pairwise,
  #     1 = Query-anchored showing identities,
  #     2 = Query-anchored no identities,
  #     3 = Flat query-anchored showing identities,
  #     4 = Flat query-anchored no identities,
  #     5 = BLAST XML,
  #     6 = Tabular,
  #     7 = Tabular with comment lines,
  #     8 = Seqalign (Text ASN.1),
  #     9 = Seqalign (Binary ASN.1),
  #    10 = Comma-separated values,
  #    11 = BLAST archive (ASN.1),
  #    12 = Seqalign (JSON),
  #    13 = Multiple-file BLAST JSON,
  #    14 = Multiple-file BLAST XML2,
  #    15 = Single-file BLAST JSON,
  #    16 = Single-file BLAST XML2,
  #    18 = Organism Report
  #
  #   Options 6, 7 and 10 can be additionally configured to produce
  #   a custom format specified by space delimited format specifiers,
  #   or by a token specified by the delim keyword.
  #    E.g.: "10 delim=@ qacc sacc score".
  #   The delim keyword must appear after the numeric output format
  #   specification.
  #   The supported format specifiers are:
  #   	    qseqid means Query Seq-id
  #   	       qgi means Query GI
  #   	      qacc means Query accesion
  #   	   qaccver means Query accesion.version
  #   	      qlen means Query sequence length
  #   	    sseqid means Subject Seq-id
  #   	 sallseqid means All subject Seq-id(s), separated by a ';'
  #   	       sgi means Subject GI
  #   	    sallgi means All subject GIs
  #   	      sacc means Subject accession
  #   	   saccver means Subject accession.version
  #   	   sallacc means All subject accessions
  #   	      slen means Subject sequence length
  #   	    qstart means Start of alignment in query
  #   	      qend means End of alignment in query
  #   	    sstart means Start of alignment in subject
  #   	      send means End of alignment in subject
  #   	      qseq means Aligned part of query sequence
  #   	      sseq means Aligned part of subject sequence
  #   	    evalue means Expect value
  #   	  bitscore means Bit score
  #   	     score means Raw score
  #   	    length means Alignment length
  #   	    pident means Percentage of identical matches
  #   	    nident means Number of identical matches
  #   	  mismatch means Number of mismatches
  #   	  positive means Number of positive-scoring matches
  #   	   gapopen means Number of gap openings
  #   	      gaps means Total number of gaps
  #   	      ppos means Percentage of positive-scoring matches
  #   	    frames means Query and subject frames separated by a '/'
  #   	    qframe means Query frame
  #   	    sframe means Subject frame
  #   	      btop means Blast traceback operations (BTOP)
  #   	    staxid means Subject Taxonomy ID
  #   	  ssciname means Subject Scientific Name
  #   	  scomname means Subject Common Name
  #   	sblastname means Subject Blast Name
  #   	 sskingdom means Subject Super Kingdom
  #   	   staxids means unique Subject Taxonomy ID(s), separated by a ';'
  #   			 (in numerical order)
  #   	 sscinames means unique Subject Scientific Name(s), separated by a ';'
  #   	 scomnames means unique Subject Common Name(s), separated by a ';'
  #   	sblastnames means unique Subject Blast Name(s), separated by a ';'
  #   			 (in alphabetical order)
  #   	sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
  #   			 (in alphabetical order)
  #   	    stitle means Subject Title
  #   	salltitles means All Subject Title(s), separated by a '<>'
  #   	   sstrand means Subject Strand
  #   	     qcovs means Query Coverage Per Subject
  #   	   qcovhsp means Query Coverage Per HSP
  #   	    qcovus means Query Coverage Per Unique Subject (blastn only)
  #   When not provided, the default value is:
  #   'qaccver saccver pident length mismatch gapopen qstart qend sstart send
  #   evalue bitscore', which is equivalent to the keyword 'std'
  #   Default = `0'
  # -show_gis
  #   Show NCBI GIs in deflines?
  # -num_descriptions <Integer, >=0>
  #   Number of database sequences to show one-line descriptions for
  #   Not applicable for outfmt > 4
  #   Default = `500'
  #    * Incompatible with:  max_target_seqs
  # -num_alignments <Integer, >=0>
  #   Number of database sequences to show alignments for
  #   Default = `250'
  #    * Incompatible with:  max_target_seqs
  # -line_length <Integer, >=1>
  #   Line length for formatting alignments
  #   Not applicable for outfmt > 4
  #   Default = `60'
  # -html
  #   Produce HTML output?
  # -sorthits <Integer, (>=0 and =<4)>
  #   Sorting option for hits:
  #   alignment view options:
  #     0 = Sort by evalue,
  #     1 = Sort by bit score,
  #     2 = Sort by total score,
  #     3 = Sort by percent identity,
  #     4 = Sort by query coverage
  #   Not applicable for outfmt > 4
  # -sorthsps <Integer, (>=0 and =<4)>
  #   Sorting option for hps:
  #     0 = Sort by hsp evalue,
  #     1 = Sort by hsp score,
  #     2 = Sort by hsp query start,
  #     3 = Sort by hsp percent identity,
  #     4 = Sort by hsp subject start
  #   Not applicable for outfmt != 0
  #
  # *** Query filtering options
  # -seg <String>
  #   Filter query sequence with SEG (Format: 'yes', 'window locut hicut', or
  #   'no' to disable)
  #   Default = `12 2.2 2.5'
  # -soft_masking <Boolean>
  #   Apply filtering locations as soft masks
  #   Default = `false'
  # -lcase_masking
  #   Use lower case filtering in query and subject sequence(s)?
  #
  # *** Restrict search or results
  # -gilist <String>
  #   Restrict search of database to list of GIs
  #    * Incompatible with:  seqidlist, taxids, taxidlist, negative_gilist,
  #   negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
  #   subject_loc
  # -seqidlist <String>
  #   Restrict search of database to list of SeqIDs
  #    * Incompatible with:  gilist, taxids, taxidlist, negative_gilist,
  #   negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
  #   subject_loc
  # -negative_gilist <String>
  #   Restrict search of database to everything except the specified GIs
  #    * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
  #   negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
  #   subject_loc
  # -negative_seqidlist <String>
  #   Restrict search of database to everything except the specified SeqIDs
  #    * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
  #   negative_gilist, negative_taxids, negative_taxidlist, remote, subject,
  #   subject_loc
  # -taxids <String>
  #   Restrict search of database to include only the specified taxonomy IDs
  #   (multiple IDs delimited by ',')
  #    * Incompatible with:  gilist, seqidlist, taxidlist, negative_gilist,
  #   negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
  #   subject_loc
  # -negative_taxids <String>
  #   Restrict search of database to everything except the specified taxonomy IDs
  #   (multiple IDs delimited by ',')
  #    * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
  #   negative_gilist, negative_seqidlist, negative_taxidlist, remote, subject,
  #   subject_loc
  # -taxidlist <String>
  #   Restrict search of database to include only the specified taxonomy IDs
  #    * Incompatible with:  gilist, seqidlist, taxids, negative_gilist,
  #   negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
  #   subject_loc
  # -negative_taxidlist <String>
  #   Restrict search of database to everything except the specified taxonomy IDs
  #    * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
  #   negative_gilist, negative_seqidlist, negative_taxids, remote, subject,
  #   subject_loc
  # -ipglist <String>
  #   Restrict search of database to list of IPGs
  #    * Incompatible with:  subject, subject_loc
  # -negative_ipglist <String>
  #   Restrict search of database to everything except the specified IPGs
  #    * Incompatible with:  subject, subject_loc
  # -entrez_query <String>
  #   Restrict search with the given Entrez query
  #    * Requires:  remote
  # -db_soft_mask <String>
  #   Filtering algorithm ID to apply to the BLAST database as soft masking
  #    * Incompatible with:  db_hard_mask, subject, subject_loc
  # -db_hard_mask <String>
  #   Filtering algorithm ID to apply to the BLAST database as hard masking
  #    * Incompatible with:  db_soft_mask, subject, subject_loc
  # -qcov_hsp_perc <Real, 0..100>
  #   Percent query coverage per hsp
  # -max_hsps <Integer, >=1>
  #   Set maximum number of HSPs per subject sequence to save for each query
  # -culling_limit <Integer, >=0>
  #   If the query range of a hit is enveloped by that of at least this many
  #   higher-scoring hits, delete the hit
  #    * Incompatible with:  best_hit_overhang, best_hit_score_edge
  # -best_hit_overhang <Real, (>0 and <0.5)>
  #   Best Hit algorithm overhang value (recommended value: 0.1)
  #    * Incompatible with:  culling_limit
  # -best_hit_score_edge <Real, (>0 and <0.5)>
  #   Best Hit algorithm score edge value (recommended value: 0.1)
  #    * Incompatible with:  culling_limit
  # -subject_besthit
  #   Turn on best hit per subject sequence
  # -max_target_seqs <Integer, >=1>
  #   Maximum number of aligned sequences to keep
  #   (value of 5 or more is recommended)
  #   Default = `500'
  #    * Incompatible with:  num_descriptions, num_alignments
  #
  # *** Statistical options
  # -dbsize <Int8>
  #   Effective length of the database
  # -searchsp <Int8, >=0>
  #   Effective length of the search space
  # -sum_stats <Boolean>
  #   Use sum statistics
  #
  # *** Search strategy options
  # -import_search_strategy <File_In>
  #   Search strategy to use
  #    * Incompatible with:  export_search_strategy
  # -export_search_strategy <File_Out>
  #   File name to record the search strategy used
  #    * Incompatible with:  import_search_strategy
  #
  # *** Extension options
  # -xdrop_ungap <Real>
  #   X-dropoff value (in bits) for ungapped extensions
  # -xdrop_gap <Real>
  #   X-dropoff value (in bits) for preliminary gapped extensions
  # -xdrop_gap_final <Real>
  #   X-dropoff value (in bits) for final gapped alignment
  # -window_size <Integer, >=0>
  #   Multiple hits window size, use 0 to specify 1-hit algorithm
  # -ungapped
  #   Perform ungapped alignment only?
  #
  # *** Miscellaneous options
  # -parse_deflines
  #   Should the query and subject defline(s) be parsed?
  # -num_threads <Integer, >=1>
  #   Number of threads (CPUs) to use in the BLAST search
  #   Default = `1'
  #    * Incompatible with:  remote
  # -mt_mode <Integer, (>=0 and =<1)>
  #   Multi-thread mode to use in BLAST search:
  #    0 (auto) split by database
  #    1 split by queries
  #   Default = `0'
  #    * Requires:  num_threads
  # -remote
  #   Execute search remotely?
  #    * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
  #   negative_gilist, negative_seqidlist, negative_taxids, negative_taxidlist,
  #   subject_loc, num_threads
  # -use_sw_tback
  #   Compute locally optimal Smith-Waterman alignments?



