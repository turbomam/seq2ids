create view smin_smax as
select
	br.*,
	min(cast(sstart as int), cast(send as int)) as smin,
	max(cast(sstart as int), cast(send as int)) as smax
from
	blast_results br ;
