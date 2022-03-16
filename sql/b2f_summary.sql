create view b2f_summary as
select
	qacc,
	sacc ,
	bitscore ,
	evalue ,
	"length" ,
	qcovs ,
	salltitles ,
	sstart ,
	send ,
	smin ,
	smax ,
	sskingdoms ,
	scomnames ,
	staxids,
	record_annotations_data_file_division,
	record_annotations_organism,
	feature_type ,
	feature_location_start ,
	feature_location_end,
	feature_qualifier_allele, feature_qualifier_bound_moiety, feature_qualifier_chromosome, feature_qualifier_clone, feature_qualifier_db_xref, feature_qualifier_EC_number, feature_qualifier_gene, feature_qualifier_inference, feature_qualifier_locus_tag, feature_qualifier_mobile_element_type, feature_qualifier_note, feature_qualifier_plasmid, feature_qualifier_product, feature_qualifier_protein_id, feature_qualifier_regulatory_class, feature_qualifier_standard_name
from
	smin_smax ss
join features f on
	ss.sacc = f.record_name
	and f.feature_location_start >= ss.smin
	and f.feature_location_end <= ss.smax
join parts_sequences_plus psp on
	ss.qacc = psp.seq_name
--where
--	psp."type" = 'insertion'
order by
	qacc ,
	cast(bitscore as int) desc,
	sacc ,
	feature_location_start ,
	feature_location_end ;
