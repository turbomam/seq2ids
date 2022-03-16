--CREATE INDEX blast_results_idx ON blast_results(qacc, sacc, bitscore, "length", pident, qcovs, sstart, send);

CREATE INDEX parts_sequences_plus_idx ON parts_sequences_plus(seq_name, "type");
