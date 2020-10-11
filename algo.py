import math

code = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
ind = {'A':0,'C':1,'G':2,'T':3}
sym = {0:'A',1:"C",2:'G',3:'T'}

# count the number of repetiotions of the given kmer in a text
def patternCount(text, pat):
	count = 0
	for i in range (0, len(text)-len(pat)+1):
		if text[i:i+len(pat)] == pat :
			count += 1
	print(count)
# print(patternCount('GCGCG', 'GCG'))

# identify the most frequent kmers in a text
def freqPat(text, k):
	kmers = {}
	max = 0
	for i in range (0, len(text)-k+1):
		if text[i:i+k] in kmers:
			kmers[text[i:i+k]] += 1
			if kmers[text[i:i+k]] > max:
				max = kmers[text[i:i+k]]
		else: 
			kmers[text[i:i+k]] = 1
	for i in kmers.keys():
		if kmers[i] == max:
			print(i)
# freqPat("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)

# reverse complimentary sequence
def reverseCompl(text):
	res = ""
	for i in text:
		res = code[i] + res
	return(res)
# reverseCompl("AAAACCCGGT")

# find indexes of repeating kmer in a text
def patternInd(pat, text):
	for i in range (0, len(text)-len(pat)+1):
		if text[i:i+len(pat)] == pat :
			print(i)
# patternInd("ATAT", "GATATATGCATATACTT")

# g - c
def skewGenome(text):
	count = 0
	for i in text:
		if i == "G":
			count += 1
		if i == "C":
			count -= 1
		print(count)
# skewGenome("CATGGGCATCGGCCATACGCC")		

# minimizing skew
def minSkew(text):
	count = 0
	minim = 0
	indexes = []
	for i in range(0, len(text)):
		if text[i] == "G":
			count += 1
		if text[i] == "C":
			count -= 1
		if count < minim:
			minim = count
			indexes.clear()
		if count == minim:
			indexes.append(i+1)  
	return(indexes)
# minSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")

def ApproximatePatternMatching(pat, text, d):
    count = []
    for i in range(0, len(text)-len(pat)+1):
        if HammingDistance(pat, text[i:i+len(pat)]) <= d:
            print(i)
            count.append(i)
    return count

def HammingDistance(seq1, seq2):
	hd = 0
	for i in range(0, len(seq1)):
		if seq1[i] != seq2[i]:
			hd += 1
	return hd

def Neighbors(Pattern, d):
	if d == 0:
		return [Pattern]
	if len(Pattern) == 1:
		return ["A","C","G","T"]
	Neighborhood = []
	SuffixNeighbors = Neighbors(Pattern[1:],d)
	for text in SuffixNeighbors:
		if HammingDistance(Pattern[1:], text) < d:
			for x in ["A","C","G","T"]:
				Neighborhood.append(x + text)
		else:
			Neighborhood.append(Pattern[0] + text)
	return Neighborhood

def PatternToNumber(kmer):
	if len(kmer) == 1:
		return ind[kmer]
	return 4*PatternToNumber(kmer[:-1]) + PatternToNumber(kmer[-1])

def NumberToPattern(ind,k):
	if k == 1:
		return sym[ind]
	prefixIndex = ind//4
	r = ind%4
	symbol = sym[r]
	prefixPattern = NumberToPattern(prefixIndex,k-1)
	return prefixPattern + symbol

def ComputingFrequenciesWithMismatches(Text, k, d):
	frequencyArray = []
	for i in range(0, 4**k):
		frequencyArray.append(0)
	for i in range(0, len(Text)-k+1):
		pattern = Text[i:i+k]
		neighborhood = Neighbors(pattern,d)
		for approx in neighborhood:
			j = PatternToNumber(approx)
			frequencyArray[j] +=1
	return frequencyArray

def FrequentWordsWithMismatches(Text, k,d):
	res = []
	frequencyArray = ComputingFrequenciesWithMismatches(Text,k,d)
	most = max(frequencyArray)
	for i in range(0, len(frequencyArray)):
		if frequencyArray[i] == most:
			res.append(NumberToPattern(i,k))
	return res

def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
	res = []
	frequencyArray = ComputingFrequenciesWithMismatchesAndReverseComplements(Text,k,d)
	most = max(frequencyArray)
	for i in range(0, len(frequencyArray)):
		if frequencyArray[i] == most:
			res.append(NumberToPattern(i,k))
	return res

def ComputingFrequenciesWithMismatchesAndReverseComplements(Text, k, d):
	frequencyArray = []
	for i in range(0, 4**k):
		frequencyArray.append(0)
	for i in range(0, len(Text)-k+1):
		pattern = Text[i:i+k]
		neighborhood = Neighbors(pattern,d)
		for approx in neighborhood:
			j = PatternToNumber(approx)
			frequencyArray[j] +=1
		reversedneighborhood = Neighbors(reverseCompl(pattern),d)
		for approx in reversedneighborhood:
			j = PatternToNumber(approx)
			frequencyArray[j] +=1
	return frequencyArray

def MotifEnumeration(Dna, k,d):
	patterns = []
	motifs = []
	for i in range(0, 4**k):
		motifs.append(0)
	for i in range(0,len(Dna[0])-k+1):
		kmer = Dna[0][i:i+k]
		neighbours =  Neighbors(kmer,d)
		for pat in neighbours:
			motifs[PatternToNumber(pat)] = 1
	for i in range(0,len(motifs)):
		if motifs[i] == 1:
			motif = NumberToPattern(i,k)
			for str in Dna[1:]:
				for j in range(0, len(str)-k+1):
					if HammingDistance(motif, str[j:j+k]) <= d:
						motifs[i] += 1
						if motifs[i] == len(Dna):
							patterns.append(motif)
						break
	return patterns

def MedianString(Dna, k):
	distance = len(Dna)+1
	median = []
	for i in range(0,4**k):
		pattern = NumberToPattern(i,k)
		score = d(pattern,Dna)
		if distance > score:
			distance = score
			median = []
			median.append(pattern)
		if distance == score:
			median.append(pattern)
	return median


def d(pattern, Dna):
	res = 0
	for text in Dna:
		distance = len(pattern)
		for i in range(0,len(text)-len(pattern)+1):
			hd = HammingDistance(pattern,text[i:i+len(pattern)])
			if hd < distance:
				distance = hd
		res += distance
	return res

def ProfileMostProbableKmer(text,k,profile):
	print(profile)
	probability = 0
	res = ''
	for i in range(0,len(text)-k+1):
		pattern = text[i:i+k]
		score = 1
		for n in range(0,len(pattern)):
			score *= profile[n][ind[pattern[n]]]
		if probability < score:
			probability = score
			res = pattern
	if res == '':
		return text[0:k]
	return res

def GreedyMotifSearch(Dna, k, t):
	bestmotifs = takefirstkmers(Dna,k,t)
	motifs = []
	for i in range(0,len(Dna[0])-k+1):
		motifs.append(Dna[0][i:i+k])
		for line in range(1,t):
			profile = formProfile(motifs,k)
			motifs.append(ProfileMostProbableKmer(Dna[line],k,profile))
		if Score(motifs) < Score(bestmotifs):
			bestmotifs = motifs
		motifs = []
	return bestmotifs

def takefirstkmers(Dna,k,t):
	matrix = []
	for i in range(0,t):
		matrix.append(Dna[i][0:k])
	return matrix

def formProfile(motifs,k):
	profile = []
	for i in range(0,k):
		profile.append([0,0,0,0])
		for motif in motifs:
			profile[i][ind[motif[i]]] += 1/len(motifs)
	return profile

def Score(matrix):
	profile = formProfile(matrix,len(matrix[0]))
	motif = ''
	for n in profile:
		motif += sym[n.index(max(n))]
	score = 0
	for kmer in matrix:
		score += HammingDistance(motif,kmer)
	# print(score)
	return score
