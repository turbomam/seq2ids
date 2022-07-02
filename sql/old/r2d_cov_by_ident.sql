create view ranges_to_download as select
	distinct
	sacc,
	NULL as minpos,
	NULL as maxpos
from
	blast_results
where
	qcovs > 90 and pident > 90
	and blast_db = 'swissprot';
