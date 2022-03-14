create view insertion_genome_list as select
	distinct sacc
	--	, count(1) as query_count
from
	insertion_querys_genome_ranking
where
	multi_rank = 1
	--group by
	--	sacc
	--order by
	--	count(1) desc
order by
	sacc
	;