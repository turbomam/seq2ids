.mode tabs
.headers on
select
-- psp."id" as seq_id,
	-- needs seq len. note all flanks are a fixed len.
	-- also needs a locus/feature identifier
	-- possibly other keywords like compound, CDSpartial
	-- also join modifications?
	-- get GO terms
	-- only getting partial parts because of columns containing carriage returns
	seq_name ,
	alias || '_' || seq_len || '_bp_' || psp."type" as seq_name_recon,
--	part_id ,
	seq_len,
	psp."type" as seq_type,
	pp."type" as mod_type,
	pp.alias,
	br.sacc,
	br.blast_db,
	br.bitscore ,
	ua."Gene names  (primary )" as gene_name_1ry,
	ua."Gene names" as gene_names,
	ua.Organism as organism,
	m.category ,
	ua."Protein names" as prot_names,
	ua."Function [CC]" as prot_function,
	ua."Gene ontology (GO)" as all_go
from
	parts_sequences_plus psp
join parts_partial pp on
	psp.part_id = pp.id
left outer join modifications m on
    pp.alias = m.element_id
left outer join blast_results br on
	psp.id = br.qacc
left OUTER join uniprot_annotations ua on
	br.sacc = ua.Entry
order by
	seq_name asc,
	bitscore desc;
