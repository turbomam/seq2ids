create view ranges_to_download as select distinct
	br.sacc,
	min(cast(sstart as int), cast(send as int)) as minpos,
	max(cast(sstart as int), cast(send as int)) as maxpos
from
	insertion_queries_best_genomes iqbg
join blast_results br on
	iqbg.sacc = br.sacc
group by
	iqbg.sacc
	order by br.sacc ;