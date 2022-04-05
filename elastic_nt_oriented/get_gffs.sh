while read genome; do
  echo "$genome"
  wget -O gff/$genome.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=$genome"
done <target/insertion_genome_list.tsv