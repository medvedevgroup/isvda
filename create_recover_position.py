"""
    input: vcf file, 1 based
    output: recover file, 0 based

    Give snp information in vcf file format(donor compared to reference),
    create recover file so that all the coordinates in donor will be recovered to coordinates in ref
    
    Work with script recovering_snp_position.py
"""
from sys import argv

if len(argv) < 2:
    print "Error: input vcf file required!\n"
    print "usage: " + argv[0] +" vcfFile(high confident SNPs) > recover_to_ref_file"
    exit()

vcfFilename = argv[1]

print "#POS\tEVENT\tLEN"
print "#POS, position before event happen, 0 based"
print "#EVENT, -: after POS, section been deleted; +: after POS, section been inserted"
print "#LEN, length of section"

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
id_pos = {}
id_event = {}
id_len = {}
idNum = 0
for pos in posList:
    refLen = pos_refLen[pos]
    altLen = pos_altLen[pos]
    previousPosition = -1
    if (refLen > altLen):
        #del happen
        previousPosition = pos + altLen - 1
        eventLen = refLen-altLen
        for storeId in id_pos:
            id_pos[storeId] -= eventLen
        id_pos[idNum] = previousPosition
        id_event[idNum] = "-"
        id_len[idNum] = eventLen
    elif (refLen < altLen):
        #ins happen
        previousPosition = pos + refLen - 1
        eventLen = altLen - refLen
        for storeId in id_pos:
            id_pos[storeId] += eventLen
        id_pos[idNum] = previousPosition
        id_event[idNum] = "+"
        id_len[idNum] = eventLen
    idNum += 1

"""
    output snp position in donor genome, in ordered sequence
""" 
idList = id_pos.keys()
idList.sort(reverse=True)
for idNum in idList:
    #output
    pos = str(id_pos[idNum])
    event = id_event[idNum]
    eventLen = str(id_len[idNum])
    print pos + "\t" + event + "\t" + eventLen
