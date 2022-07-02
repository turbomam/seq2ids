import re

import requests
from bs4 import BeautifulSoup

url = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
txt_out_fp = '../nt_file_list.txt'

req = requests.get(url)
soup = BeautifulSoup(req.text, "html.parser")

file_list = []

for i in soup.findAll('a'):
    file_list.append(i.text)

r = re.compile("^nt")

filtered_list = list(filter(r.match, file_list))

# todo with

textfile = open(txt_out_fp, "w")

for element in filtered_list:
    textfile.write(element + "\n")

textfile.close()
