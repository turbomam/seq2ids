attach 'local/live_sqlite.db' as felix;
attach 'target/seq2ids.db' as seq2ids;
DROP TABLE IF EXISTS seq2ids.parts_sequences_plus;
create table seq2ids.parts_sequences_plus as select * from felix.parts_sequences;
create table seq2ids.parts_partial as select * from felix.parts;
create table seq2ids.modifications as select * from felix.modifications;
update seq2ids.parts_sequences_plus set "sequence" = replace(replace("sequence", ' ', ''), '\n', '');
ALTER TABLE seq2ids.parts_sequences_plus ADD seq_len int;
update seq2ids.parts_sequences_plus set seq_len = length("sequence");
