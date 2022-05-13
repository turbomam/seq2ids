create view parts_up_annotations as
select
	psp.id as seq_id,
	psp.part_id ,
	psp.seq_name ,
	psp."type" ,
	cast(psp.seq_len as int ) as seq_len_int,
	br.evalue ,
	br.bitscore ,
	br.qcovs ,
	br."length" * 3 as hsp_length_int_x3,
	br.pident ,
	br.sacc,
	ua."Entry name" as entry_name,
	ua."Gene names" as gene_names,
	ua."Gene names  (primary )" as primary_names ,
	ua."Gene names  (synonym )" as synonyms,
	ua."Gene names  (ORF )" as orf_names,
	ua.Organism ,
	ua."Protein names" as protein_names,
	ua."Cross-reference (BRENDA)" as brenda,
	ua."Function [CC]" as up_function,
	ua."EC number" as ec,
	ua.Annotation,
	ua."Biotechnological use" as up_biotech,
	ua."Gene ontology (biological process)" as bp,
	ua."Gene ontology (molecular function)" as mf
from
	blast_results br
join uniprot_annotations ua on
	sacc = ua.Entry
-- join parts_sequences_plus psp on
--	qacc = psp.seq_name
-- todo double CHECK
-- todo cast psp.id to string?
join parts_sequences_plus psp on
	br.qacc = psp.id
order by seq_name, bitscore desc;
