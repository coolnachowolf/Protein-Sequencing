"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file = open(filename, "r")
    words = ''
    for line in file:
        if len(line) > 1:
            line = line.strip()
            words += line
    file.close()
    return words


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    dna_lst = []
    s = ''
    rna_lst = []
    end_value = ["TAA", "TAG", "TGA"]
    i = startIndex
    while i in range(len(dna)):
        s += dna[i]
        i += 1
        if(len(s)==3):
            dna_lst.append(s)
            if(s in end_value):
                break
            else:
                s = ''  
    for string in dna_lst:
        rna_lst.append(string.replace("T", "U"))
    return rna_lst


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f = open(filename)
    data = json.load(f)
    codon_dict = {}
    for amino_value in data:
        for value in data[amino_value]:
            codon_value = ''
            for i in range(len(value)):
                if(value[i]=="T"):
                    codon_value += "U"
                else:
                    codon_value += value[i]
            codon_dict[codon_value] = amino_value
    return codon_dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein_lst = []
    end_lst = ["UAA", "UAG", "UGA"]
    for i in range(len(codons)):
        if codons[i] == "AUG" and len(protein_lst)==0:
            protein_lst.append("Start")
            if codons[i] in end_lst:
                protein_lst.append("Stop")
        else:
            protein_lst.append(codonD[codons[i]])
    # print("protein_lst=", protein_lst)
    return protein_lst


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna = readFile(dnaFilename)
    codonD = makeCodonDictionary(codonFilename)
    proteins_lst = []
    unused_count = 0
    loop = 0
    while loop in range(len(dna)):
        if dna[loop : loop+3]=="ATG":
            codons = dnaToRna(dna,loop)
            protein = generateProtein(codons,codonD)
            proteins_lst.append(protein)
            loop += 3*len(protein)
        else:
            loop += 1
            unused_count += 1
    print(len(proteins_lst),unused_count)
    return proteins_lst


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    unique_lst = []
    for i in proteinList1:
        if i not in unique_lst:
            unique_lst.append(i)
    common_lst = []
    for i in unique_lst:
        if i in proteinList2:
            if i not in common_lst:
                common_lst.append(i)
    return common_lst


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combine_lst = []
    for i in range(len(proteinList)):
        for j in range(len(proteinList[i])):
            combine_lst.append(proteinList[i][j])
    return combine_lst


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    return


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    return


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":

    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()
    

    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()


    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
