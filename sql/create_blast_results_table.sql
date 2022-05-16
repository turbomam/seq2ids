-- blast_results definition

CREATE TABLE "blast_results"(
  "qacc" TEXT,
  "qcovhsp" REAL,
  "qcovs" REAL,
  "qstart" INTEGER,
  "qend" INTEGER,
  "bitscore" REAL,
  "score" REAL,
  "evalue" REAL,
  "length" INTEGER,
  "pident" REAL,
  "sacc" TEXT,
  "sstart" INTEGER,
  "send" INTEGER,
  "sallacc" TEXT,
  "sseqid" TEXT,
  "sallseqid" TEXT,
  "stitle" TEXT,
  "salltitles" TEXT,
  "staxids" TEXT,
  "sskingdoms" TEXT,
  "sscinames" TEXT,
  "sblastnames" TEXT,
  "scomnames" TEXT,
  "blast_db" TEXT
);

CREATE INDEX blast_results_idx ON blast_results(qacc, sacc, bitscore, "length", pident, qcovs, sstart, send);
