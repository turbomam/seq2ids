create view insertion_genome_hsp_counts as
select
	sacc,
	count(1) as hsp_count
from
	blast_results br
join parts_sequences_plus psp
	on
	br.qacc = psp.seq_name
--where
--	psp."type" = 'insertion'
group by
	sacc
order by
	count(1) desc;


