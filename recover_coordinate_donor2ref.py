"""
    input:
    1.vcf file describe snp info
    2.recover file that has been generated, 0-based
    output: bed file about snp position and length, 0-based

    Give vcf file that contains coordinates in donor sequence(in this script, it contains snp coordinate info)
    and the recover file of reference(containing coordinates changing information, this can be acquired by script create_recover_position.py)

    Generate corresponding coordinates in reference genome.

    Pay attention that the corresponding coordinate may not exist on reference genome because insertion happened on reference.
    Here we did not deal with this situation.

    Author: Chen Sun(cxs1031@cse.psu.edu)
    Modified by Chen Sun on 11/6/2014, to correctly recover.
    The idea is to start sorted, not reverse sorted.

    Modified by Chen Sun on 12/10/2014, to output vcf file format.

    In the future, we need to modified this so that it can modify a region, not just a point, will contain start and end.
"""
import sys,getopt

def find_position(recoverPos, snpPosList):
    startP = 0
    endP = len(snpPosList) - 1
    midP = int((startP + endP)/2)
    while(endP-startP > 1):
        if(recoverPos < snpPosList[midP]):
            endP = midP
            midP = int((startP + endP)/2)
        elif(recoverPos > snpPosList[midP]):
            startP = midP
            midP = int((startP + endP)/2)
        else:
            return midP
    # there may not be exact position, we just want to find a search start point
    # so to avoid error, we start from smaller position
    return startP
    
def usage():
    print "Usage: " + sys.argv[0] + "\n\t[-v vcf_file_from_donor | -b bed_file_from_donor] \n\t-r recover_file_from_reference \n\t[-f output_file_in_vcf | -d output_file_in_bed]"

def main(argv):
    if len(argv) < 1:
        usage()
        sys.exit()
    
    vcfFilename = ""
    fileFormat = "vcf" # it may indicate it is vcf file, but if bed file, I will indicate it.
    recoverFilename = ""
    outputFormat = "vcf"
    outputFilename = ""
    
    try:
        opts, args = getopt.getopt(argv, "hv:b:r:f:d:", [])
    except getopt.GetoptError as err:
        print err
        sys.exit()

    for opt, arg in opts:
        if opt in ("-v"):
            vcfFilename = arg
        if opt == "-h":
            usage()
            sys.exit()
        if opt in ("-b"):
            vcfFilename = arg
            fileFormat = "bed"
        if opt in ("-r"):
            recoverFilename = arg
        if opt in ("-f"):
            outputFilename = arg
        if opt in ("-d"):
            outputFilename = arg
            outputFormat = "bed"
    
    snpPosList = []
    newSnpList = []
    
    vcfFile = open(vcfFilename)
    
    headInfo = ""
    pos_originalInfo = {}
    for line in vcfFile.readlines():
        if line.startswith("#") :
            headInfo += line;
            continue
        
        line = line.strip()
        columns = line.split("\t")
        pos = int(columns[1])-1 # here we trans 1-based to 0 based 
        if fileFormat == "bed": # by doing this, we can deal with both vcf or bed
            pos += 1
        snpPosList.append(pos)
        pos_originalInfo[pos] = line
    vcfFile.close() #==========close vcf file
    snpPosList.sort()
    newSnpList = snpPosList[:] #copy snp list to a new list by value
    
    
    recoverFile = open(recoverFilename)
    
    pos_event = {}
    pos_len = {}
    for line in recoverFile.readlines():
        if line.startswith("#"):
            continue
        line = line.strip()
        columns = line.split("\t")
        pos = int(columns[0])
        event = columns[1]
        eventLen = int(columns[2])
        pos_event[pos] = event
        pos_len[pos] = eventLen
    recoverFile.close() #======close recover file

    recoverPosList = pos_event.keys()
    recoverPosList.sort(reverse=True)

    for recoverPos in recoverPosList:
        event = pos_event[recoverPos]
        eventLen = pos_len[recoverPos]
        change = 0
        if(event == "+"):
            change = eventLen * -1
        elif(event == "-"):
            change = eventLen
        
        #binary search in snpPosList
        index = find_position(recoverPos, snpPosList)
        for i in range(index, len(snpPosList)):
            if recoverPos < snpPosList[i]:
                for j in range(i, len(snpPosList)):
                    newSnpList[j] += change #snpPosList is just for comparing, newSnpList store the real value
                break
            
    # newSnpList.sort()

    outputFile = open(outputFilename,"w")

    #output a bed format mark
    if outputFormat == "bed":
        newSnpList.sort()
        outputFile.write("#bed\n")
        for position in newSnpList:
            outputFile.write("chr17\t" + str(position) + "\t" + str(position+1)+"\n")
    elif outputFormat == "vcf":
        outputFile.write(headInfo)
        for i in range(len(snpPosList)):
            oldPos = snpPosList[i]
            newPos = newSnpList[i] + 1 
            originalInfo = pos_originalInfo[oldPos]
            columns = originalInfo.split("\t")
            columns[1] = str(newPos)
            newInfo = "\t".join(columns)
            newInfo += "\n"
            outputFile.write(newInfo)

    outputFile.close()

if __name__ == '__main__':
    main(sys.argv[1:])
