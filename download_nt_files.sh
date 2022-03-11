while read line; do
  echo "$line"
  ~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -T -l300m  anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/$line .
done <nt_file_list.txt
