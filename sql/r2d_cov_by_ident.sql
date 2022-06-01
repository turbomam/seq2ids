create view ranges_to_download as select
	distinct
	sacc,
	NULL as minpos,
	NULL as maxpos
from
	blast_results
where
    -- only get annotations for blast results if the qcovs and pident are both roughly > 90
	qcovs * pident > 8100
	and blast_db = 'swissprot';
