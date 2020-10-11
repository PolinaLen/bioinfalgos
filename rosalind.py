ind = {'A':0,'C':1,'G':2,'T':3}
def countingNucleotides(Dna):
    count = [0,0,0,0]
    for n in Dna:
        count[ind[n]] +=1
    return count

def transcribeToRna(Dna):
    Rna = ''
    for n in Dna:
        if n == 'T':
            Rna += 'U'
        else: 
            Rna += n
    return Rna