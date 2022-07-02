# --    - SC_curated_enzyme_name
# --    - SC_curated_gene_name
# --    - SC_curated_protein_name
# --    - SC_curated_uniprot_id
# --    - seq2ids_gene_name_1ry
# --    - seq2ids_prot_function
# --    - seq2ids_prot_names
import sqlite3

import pandas as pd

seq2ids_db_fp = "target/seq2ids.db"

worthiness_tsv_out = "target/worthiness.tsv"

worthiness_q = """
select
    "alias" as mod_alias,
    "descriptor",
    "element",
    "name",
    "notes",
    "position",
    "size",
    "status",
    aa_change,
    annotationScore * pident * qcovs / 10000 as worthiness,
    annotationScore,
    bio_safety_level,
    bitscore,
    blast_db,
    br."length" as alingment_length,
    category,
    creator_id,
    ecNumbers,
    element_id,
    element_id_number,
    element_name,
    element_species,
    geneName,
    go_ids,
    go_labels,
    modification_type,
    org_scientificName,
    pident,
    primaryAccession,
    principal_investigator,
    principal_investigator_email as email,
    psp."type" as seq_type,
    psp.id as seq_id,
    psp.seq_len,
    psp.seq_name,
    qcovs,
    rec_full_name,
    sacc,
    subcategory_size,
    taxon_id,
    taxonId,
    uniProtkbId
from
    modifications m
join parts_partial pp 
on
    m.parts_ptr_id = pp.id
join parts_sequences_plus psp 
on
    pp.id = psp.part_id
join blast_results br 
on
    psp.id = br.qacc
left join uniprot_annotations ua 
on
    br.sacc = ua.primaryAccession ;
"""

seq2ids_db_conn = sqlite3.connect(seq2ids_db_fp)

worthiness_res = pd.read_sql_query(worthiness_q, seq2ids_db_conn)

worthiness_res.loc[worthiness_res['blast_db'].eq("fpbase"), "worthiness"] = worthiness_res.loc[
                                                                                worthiness_res['blast_db'].eq(
                                                                                    "fpbase"), "pident"] * \
                                                                            worthiness_res.loc[
                                                                                worthiness_res['blast_db'].eq(
                                                                                    "fpbase"), "qcovs"]

# todo get more FPbase details
# todo won't report all insertions for compound insertions
worthiness_res["rank"] = worthiness_res.groupby(["seq_name"])["worthiness"].rank("dense", ascending=False)

worthiest = worthiness_res.loc[worthiness_res['rank'].eq(1)]

worthiest.to_csv(worthiness_tsv_out, sep="\t", index=False)

worthiest.to_sql(
    name="worthiest", con=seq2ids_db_conn, if_exists="replace", index=False
)

# --funding_source
# --genotype_phenotype
# --keywords
# --summary
# --references
# --intellectual_property

#
