#   aminocount.py
#   Input: FASTA file containing translated CDS sequences
#   Output: total count of amino acids in the entire file and a CSV file of the count for individual CDS

import csv
from collections import Counter
from Bio import SeqIO

readFile = open("RAP-DB_protein.fasta", 'r')
writeFileName = "RAP-DB_protein.csv"

data = ""
header = ["GenID", "Length"]
records = []

# count total number of amino acids to arrange csv in descending order of abundance
for line in readFile:
	if line[0] != ">":
		data = data + line
data = data.replace("\n", "").replace("*", "")
aminoCount = Counter(data)
for item in aminoCount.most_common():
	header.append(item[0])

# reset file pointer and count for each sequence
readFile.seek(0)
sequences = SeqIO.parse(readFile, "fasta")
for sequence in sequences:
	record = [sequence.id, len(sequence.seq) - 1]
	recordAminoCount = Counter(sequence.seq)
	for item in aminoCount.most_common():
		record.append(recordAminoCount[item[0]])
	records.append(record)
readFile.close()

with open(writeFileName, "w", newline="") as writeFile:
	writer = csv.writer(writeFile)
	writer.writerow(header)
	writer.writerows(records)

print(aminoCount.most_common())
