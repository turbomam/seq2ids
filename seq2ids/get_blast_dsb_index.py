from bs4 import BeautifulSoup

import requests

import re

url = "https://ftp.ncbi.nlm.nih.gov/blast/db/"

req = requests.get(url)
soup = BeautifulSoup(req.text, "html.parser")

file_list = []

for i in soup.findAll('a'):
    file_list.append(i.text)

r = re.compile("^nt")

filtered_list = list(filter(r.match, file_list))

# print(filtered_list)

textfile = open("../nt_file_list.txt", "w")

for element in filtered_list:
    textfile.write(element + "\n")

textfile.close()
