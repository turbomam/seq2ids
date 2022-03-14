create view insertion_querys_genome_ranking as
SELECT
	qacc,
	br.sacc,
	Row_number() OVER (PARTITION BY qacc
ORDER BY
	bitscore DESC,
	"length" desc,
	pident desc,
	qcovs desc,
	hsp_count desc) AS multi_rank
FROM
	blast_results br
join insertion_genome_hsp_counts ighc on
	br.sacc = ighc.sacc;
