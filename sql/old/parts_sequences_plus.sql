select
	"id" ,
	part_id ,
	"file" ,
	seq_name ,
	"type" ,
	date_added ,
	"sequence" ,
	length("sequence") as seq_len
from
	parts_sequences ps ;
