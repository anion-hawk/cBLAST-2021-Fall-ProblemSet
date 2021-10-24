#   cds.py
#   Input: Genbank file of an organism
#   Output: FASTA file containing all the CDS present in the genome

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

readFile = open("U00096.3.gb", "r")
writeFile = open("U00096.3_cds.fasta", "w")
sequences = []
count = 0

for record in SeqIO.parse(readFile, "genbank"):
    for feature in record.features:
        if feature.type.casefold() == "CDS".casefold():
            count = count + 1
            descStr = "[gene=" + feature.qualifiers["gene"][0] + "] " + "[locus_tag=" + feature.qualifiers["locus_tag"][0] + "] " + "[protein=" + feature.qualifiers["product"][0] + "] "
            if feature.qualifiers.get("pseudo") is None:
                descStr = descStr + "[protein_id=" + feature.qualifiers["protein_id"][0] + "] "
            else:
                descStr = descStr + "[pseudo=true] "
            location = feature.location.__str__().replace("[", "").replace("]", "").replace(":", "...").replace("{", "(").replace("}", ")")
            descStr = descStr + "[location=" + location + "]"
            sequence = SeqRecord(feature.extract(record.seq),
                                 id="lcl|U00096.3_cds_" + str(count),
                                 description=descStr)
            sequences.append(sequence)
print("Number of Coding Regions: " + str(count))
SeqIO.write(sequences, writeFile, "fasta")
readFile.close()
writeFile.close()
