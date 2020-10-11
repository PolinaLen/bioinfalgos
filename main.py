import algo

# file = open("dataset.txt","r")
# k,t = file.readline()[:-1].split(' ')
# Dna = []
# for line in file.readlines():
#     Dna.append(line[:-1])
# file.close() 
# # print(Dna)

Dna = ["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
"GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
"GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"]
motifs = algo.MedianString(Dna,7)
print(motifs)
