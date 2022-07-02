create view insertion_queries_best_genomes as select
	qacc ,
	sacc ,
	multi_rank
from
	(
	SELECT
		qacc,
		br.sacc,
		Row_number() OVER (PARTITION BY qacc
	ORDER BY
		cast(bitscore as int) DESC,
		--	cast(qcovs as int) desc,
		--	cast("length" as int) desc,
		--	cast(pident as int) desc,
		hsp_count desc) AS multi_rank
	FROM
		blast_results br
	join insertion_genome_hsp_counts ighc on
		br.sacc = ighc.sacc
--	join parts_sequences_plus psp on
--		br.qacc = psp.seq_name
--	where
--		"type" = 'insertion'
)
where
	multi_rank <= 10 ;