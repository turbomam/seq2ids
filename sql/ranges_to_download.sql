create view ranges_to_download as select
	br.sacc,
	min(cast(sstart as int), cast(send as int)) as minpos,
	max(cast(sstart as int), cast(send as int)) as maxpos
from
	insertion_genome_list igl
join blast_results br on
	igl.sacc = br.sacc
group by
	igl.sacc;