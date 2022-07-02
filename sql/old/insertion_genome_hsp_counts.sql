create view insertion_genome_hsp_counts as
select
	sacc,
	count(1) as hsp_count
from
	blast_results br
join parts_sequences_plus psp
	on
	br.qacc = psp.id
where br.blast_db = 'swissprot'
--	psp."type" = 'insertion'
group by
	sacc
order by
	count(1) desc;


