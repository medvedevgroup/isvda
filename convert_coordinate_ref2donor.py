"""
    input: vcf file, 1 based
           bed file, 0 based
    output: recover file, 0 based

    Given snp information that stored in vcf file, the snp is donor compared to reference, and coordinate is in reference
    bed file containing coordinate info in reference

    Convert coordinate to corresponding one on donor genome

    Pay attention that the coordinate may not exist on donor genome because of deletion happened on reference,
    here we did not deal with it.
"""
from sys import argv

if len(argv) < 3:
    print "Error: parameter required!\n"
    print "usage: " + argv[0] +" vcfFile bedFile > new_bedFile"
    exit()

vcfFilename = argv[1]
bedFilename = argv[2]

vcfFile = open(vcfFilename)
pos_refLen = {}
pos_altLen = {}

# readfile
for line in vcfFile.readlines():
    if(line.startswith("#")):
        continue
    line = line.strip()
    columns = line.split("\t")
    pos = int(columns[1]) - 1 #vcf file is 1 based
    ref = columns[3]
    alt = columns[4]
    refLen = len(ref)
    altLen = len(alt)
    
    if (refLen == altLen):
        continue
    pos_refLen[pos] = refLen
    pos_altLen[pos] = altLen

vcfFile.close()

"""
    Traverse all snp position, start from largest position
    For one position, change the coordinate after it.
"""
posList = pos_refLen.keys()
posList.sort(reverse=True)

"""
    Read bed file position.
    
"""
bedFile = open(bedFilename)

bedId = 0
bedId_bedStart = {}
bedId_bedEnd = {}
for line in bedFile.readlines():
    bedId += 1
    if line.startswith("#"):
        continue
    line = line.strip()
    columns = line.split("\t")
    startPosition = int(columns[1])
    endPosition = int(columns[2])
    bedId_bedStart[bedId] = startPosition
    bedId_bedEnd[bedId] = endPosition
bedFile.close()


for bedId in bedId_bedStart:
    startPosition = bedId_bedStart[bedId]
    for pos in posList:
        if pos > startPosition:
            continue
        refLen = pos_refLen[pos]
        altLen = pos_altLen[pos]
        if (refLen > altLen):
            #del happen
            eventLen = refLen-altLen
            bedId_bedStart[bedId] -= eventLen
        elif (refLen < altLen):
            #ins happen
            eventLen = altLen - refLen
            bedId_bedStart[bedId] += eventLen
    endPosition = bedId_bedEnd[bedId]        
    for pos in posList:
        if pos > endPosition:
            continue
        refLen = pos_refLen[pos]
        altLen = pos_altLen[pos]
        if (refLen > altLen):
            #del happen
            eventLen = refLen-altLen
            bedId_bedEnd[bedId] -= eventLen
        elif (refLen < altLen):
            #ins happen
            eventLen = altLen - refLen
            bedId_bedEnd[bedId] += eventLen
        
for bedId in bedId_bedStart:
    print "chr17\t" + str(bedId_bedStart[bedId]) + "\t" + str(bedId_bedEnd[bedId])
