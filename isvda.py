#!/usr/bin/python
"""
Author: Chen Sun, chensun@cse.psu.edu

"""
import sys
versionError = "You are using an old version of python, please upgrade to python 2.7+\n"
if sys.hexversion < 0x02070000:
    print (versionError)
#elif sys.hexversion > 0x03000000:
#    print ("python 3")

import subprocess, argparse, os
from datetime import datetime
from multiprocessing import Pool
import time

citation = "Please cite our paper."

DEBUG = False
RUN = True

def shell_run(command):
    if not RUN:
        time.sleep(5.7)
        print (command)
    else:
        print (command)
        print (".")
        subprocess.call(command, shell=True)

def serial_exe(command_list):
    #print 'running serial_exe: ', command_list
    for command in command_list:
        shell_run(command)

def check_command(command): 
    """
    check if corresponding command available
    """
    if os.path.isfile(command):
        return True

    for cmdpath in os.environ['PATH'].split(':'):
        if os.path.isdir(cmdpath) and command in os.listdir(cmdpath):
            return True
    return False

parser = argparse.ArgumentParser(epilog = citation)
parser.add_argument('-r', '--reference', required=not DEBUG, default = '/gpfs/home/cxs1031/standard/repairing/bacteria/ref.fa', help = 'reference genome file path')
parser.add_argument('-p', '--pair_end', required=not DEBUG, nargs = 2, default = '/gpfs/home/cxs1031/standard/repairing/bacteria/donor_10x_1.fq /gpfs/home/cxs1031/standard/repairing/bacteria/donor_10x_2.fq'.split(), help = 'two pair end read files')
parser.add_argument('-i', '--iteration', required=not DEBUG, default = 4, type = int, help = 'number of iteration to run, if not detect new high confident snp, the pipeline will terminate automatically, even there are still iterations left')
parser.add_argument('-w', '--workspace', required=not DEBUG, default='/gpfs/home/cxs1031/standard/repairing/bacteria', help = 'directory where result will be put in, should be empty')
parser.add_argument('-o', '--output', help = 'snp result file in vcf file format')
parser.add_argument('-t', '--thread', default=1, type=int, help = 'multi-thread number to run each step that can be paralleled')
#parser.add_argument('-k', action='store_true')
parser.add_argument('--bwa', default='bwa', help = 'bwa program path')
parser.add_argument('-e', help='use bowtie instead of bwa, default is deactivte', action='store_true')
parser.add_argument('--bowtie', default='bowtie2', help = 'bowtie2 program path')
parser.add_argument('--bindex', default='bowtie2-build', help = 'bowtie2-build program path')
parser.add_argument('--samtools', default='samtools', help = 'samtools program path')
parser.add_argument('--python', default='python', help = 'python program path')
parser.add_argument('--parallel_freebayes', required=not DEBUG, default='/gpfs/home/cxs1031/src/freebayes/scripts/', help = 'parallel freebayes directory')
parser.add_argument('--freebayes', default='freebayes', help = 'sequential freebayes program path')
parser.add_argument('-g', '--whole_genome', help='activate whole genome analysis, use when you have multi chromosome', action='store_true')
args = parser.parse_args()  

def error():
    parser.print_help()
    sys.exit()

def main():
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    print args
    
    scriptPath = sys.path[0]
    #check that there is no blank in all file path and file name
    if " " in args.reference or " " in args.pair_end:
        print ("\nError: file path and file name should not contian blank space\n")
        error()

    #check if workspace is empty
    if not os.path.exists(args.workspace):
        os.makedirs(args.workspace)

    if not check_command(args.python):
        print("\nError: python: " + args.python + " not found\n")
        error()

    if not check_command(args.samtools):
        print("\nError: samtools: " + args.samtools + " not found\n")
        error()
    
    print args.thread

    if args.thread > 1:
        #print args.thread, os.environ['MALLOC_ARENA_MAX']
        if int(args.thread) > int(os.environ['MALLOC_ARENA_MAX']):
            print ('\nError: detect only '+ os.environ['MALLOC_ARENA_MAX'] + ' threads available, but thread setting is ' + str(args.thread) + '\n')
            print ('\tPlease set thread number not greater than available number, or reset environment variable MALLOC_ARENA_MAX\n')
            error()

    #exit()

    if args.thread > 1 and not os.path.isdir(args.parallel_freebayes):
        if args.parallel_freebayes == "/gpfs/home/cxs1031/src/freebayes/scripts/":
            print ("\nError: parallel_freebayes directory not set\n")
            error()
        else:
            print ("\nError: parallel_freebayes director: " + args.parallel_freebayes + " not exist\n")
            error()

    genomeList=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    #if os.listdir(args.workspace) != []:
    #    print ("\nError: workspace should be empty\n")
    #    parser.print_help()
    #    sys.exit()

    # create fold for temp files
    for i in range(args.iteration):

        print ("running iteration " + str(i+1) + ": " + str(datetime.now()))

        iterFold = args.workspace + "/" + str(i+1)
        
        referenceFile = ""
        if i == 0:
            referenceFile = args.reference
        else:
            referenceFile = args.workspace + "/" + str(i)+"/"+str(i)+".fa"
        
        #shell_run("mkdir " + iterFold)
        
        if not os.path.exists(iterFold):
            os.makedirs(iterFold)

        #check ref.fa.fai
        refFai = referenceFile+".fai"
        if not os.path.isfile(refFai):
            faidxCommand = args.samtools + " faidx " + referenceFile
            shell_run(faidxCommand)
            #subprocess.call(faidxCommand, shell=True)

        sortPrefix = iterFold + "/" + str(i+1) + "_sort"
        
        if not args.e:
            #run bwa
            if not os.path.isfile(referenceFile+".bwt"):
                indexFastaCommand = args.bwa + " index " + referenceFile
                shell_run(indexFastaCommand)
                #subprocess.call(indexFastaCommand, shell=True)
        
            bwaThread = ""
            if args.thread > 1:
                bwaThread = " -t " + str(args.thread)
                #bwaThread = " -t 4"

            if i == 0 and not check_command(args.bwa):
                print ("\nError: bwa command: " + args.bwa + " not found. You can set by --bwa parameter.\n")
                error()
            
            bwaCommand = args.bwa + " mem" + bwaThread + " " + referenceFile + " " + args.pair_end[0] + " " + args.pair_end[1] + " | " + args.samtools + " view -Shu - | " + args.samtools + " sort - " + sortPrefix
            #bwaCommand = args.bwa + " mem" + bwaThread + " " + referenceFile + " " + args.pair_end[0] + " " + args.pair_end[1] + " | " + args.samtools + " view -Shu - > "+ iterFold + "/" + str(i+1) + ".bam"
            #| " + args.samtools + " sort - " + sortPrefix
            shell_run(bwaCommand)
            #sortCommand = args.samtools + " sort -@ " + str(args.thread) + " " + iterFold + "/" + str(i+1) + ".bam " + sortPrefix
            #shell_run(sortCommand)
        else:
            # run bowtie instead of bwa

            samFile = iterFold + "/" + str(i+1) + ".sam"
            bowtieIndexFile = referenceFile+"_bowtie2"
            bowtieIndexCommand = args.bindex + " " + referenceFile + " " + bowtieIndexFile
            if not os.path.isfile(bowtieIndexFile):
                shell_run(bowtieIndexCommand)
            bowtieCommand = args.bowtie + " -p " + str(args.thread) + " -x " + bowtieIndexFile + " -1 " + args.pair_end[0] + " -2 " + args.pair_end[1] + " -S " + samFile
            shell_run(bowtieCommand)
            #subprocess.call(bwaCommand, shell=True)

            bamFile = iterFold + "/" + str(i+1) + ".bam"
            sam2bamCommand = args.samtools + " view -Sb " + samFile + " > " + bamFile
            shell_run(sam2bamCommand)
            #subprocess.call(sam2bamCommand, shell=True)

            sortBamCommand = args.samtools + " sort " + bamFile + " " + sortPrefix
            shell_run(sortBamCommand)
            #subprocess.call(sortBamCommand, shell=True)

        sortBam = sortPrefix + ".bam"
        indexBamCommand = args.samtools + " index " + sortBam
        shell_run(indexBamCommand)
        #subprocess.call(indexBamCommand, shell=True)
        
        regionFile = args.workspace + "/" + str(i+1) + "/region_file"
        regionCommand = args.python + " " + args.parallel_freebayes + "/fasta_generate_regions.py " + refFai + " 100000 > " + regionFile
        
        vcfFile = args.workspace + "/" + str(i+1) + "/" + str(i+1) + ".vcf"
        parallelFreebayesCommand = args.parallel_freebayes + "/freebayes-parallel " + regionFile + " " + str(args.thread) + " -f " + referenceFile + " " + sortBam + " > " + vcfFile

        sequentialFreebayesCommand = args.freebayes + " -f " + referenceFile + " " + sortBam + " > " + vcfFile 
        if args.thread > 1:
            shell_run(regionCommand)
            #subprocess.call(regionCommand, shell=True)
            shell_run(parallelFreebayesCommand)
            #subprocess.call(parallelFreebayesCommand, shell=True)
        else:
            if not check_command(args.freebayes):
                print ("\nError: freebayes: " + args.freebayes + " not found\n")
                error()
            shell_run(sequentialFreebayesCommand)
            #subprocess.call(sequentialFreebayesCommand, shell=True)

        fastaFile = iterFold + "/" + str(i+1)+".fa"
        highConfidentFile = iterFold + "/" + str(i+1) + '_high_confident_on_'+str(i)+'.vcf'
        lowConfidentFile = iterFold + "/" + str(i+1) + '_low_confident_on_'+str(i)+'.vcf'

        if i == 0 and not os.path.isfile(scriptPath+"/modify_genome_with_small_variants"):
            print ("\nError: modify_genome_with_small_variants not found\n")
            error()
        if not args.whole_genome:
            constructCommand = scriptPath + "/modify_genome_with_small_variants -r " + referenceFile + " -s " + vcfFile + " -d " + fastaFile + " -h " + highConfidentFile + " -w " + lowConfidentFile
            shell_run(constructCommand)
            #subprocess.call(constructCommand, shell=True)
            recover_creating_command = args.python + ' ' + scriptPath + '/create_recover_position.py ' + iterFold + '/' + str(i+1) + '_high_confident_on_'+str(i)+'.vcf > ' + iterFold + '/' + str(i+1) + '_to_' + str(i)
            shell_run(recover_creating_command)
            
            if i > 0:
                for target in range(i,0,-1):
                    previous_iter_fold = args.workspace + "/" + str(target)
                    recovering_command = args.python + ' ' + scriptPath + '/recover_coordinate_donor2ref.py -v ' + iterFold + '/' + str(i+1) + '_high_confident_on_' + str(target) + '.vcf -r '+ previous_iter_fold + '/' + str(target) + '_to_' + str(target-1) + ' -f ' + iterFold + '/' +str(i+1) + '_high_confident_on_' + str(target-1) + '.vcf'
                    shell_run(recovering_command)

        else:
            if RUN:
                #split ref file
                genomeList = []
                REF = open(referenceFile)
                splitREF = [None for x in range(24)]
                fileIter = -1
                for line in REF.readlines():
                    if line.startswith(">"):
                        fileIter += 1
                        if fileIter > 0:
                            splitREF[fileIter-1].close()
                        line = line.strip()
                        columns = line.split(">")
                        filename = iterFold + "/ref_" + columns[1] + ".fa"
                        genomeList.append(columns[1].split('chr')[1])
                        splitREF[fileIter] = open(filename, "w")
                        splitREF[fileIter].write(">" + columns[1] + "\n")
                        continue
                    splitREF[fileIter].write(line)

                if fileIter < 23:
                    print ('\nWarning: only deal with following ' + chr(fileIter+1) + ' chromosomes\n')
                    print genomeList
                    splitREF[fileIter].close()
                    #print ('\nError: only one chromosome found in whole genome mode\n')
                    #error()
                else:
                    splitREF[fileIter].close()

                #split vcf file
                splitVCF = [None for x in range(24)]
                fileIter = 0
                chrName_iter = {}
                for g in genomeList:
                    chrName = "chr"+g
                    filename = iterFold + "/" + chrName + '_' + str(i+1) + '_on_' + str(i) + ".vcf"
                    splitVCF[fileIter] = open(filename, "w")
                    chrName_iter[chrName] = fileIter
                    fileIter += 1

                VCF = open(vcfFile)
                vcfHead = ""
                writeHead = False
                for line in VCF.readlines():
                    if line.startswith("#"):
                        vcfHead += line
                        continue
                    if not writeHead:
                        for j in range(fileIter):
                            splitVCF[j].write(vcfHead)
                        writeHead = True
                    columns = line.split("\t")
                    chrName = columns[0]
                    splitVCF[chrName_iter[chrName]].write(line)

                for j in range(fileIter):
                    splitVCF[j].close()
                
            pool = Pool(processes=args.thread)
            #construct donor genome
            print ("parallel constructing chromosomes...\n")
            for g in genomeList:
                single_command_list = []
                singleRef = iterFold + "/ref_chr" + g + ".fa"
                #singleVcf = iterFold + "/chr"+g+".vcf"
                singleVcf = iterFold + '/chr' + g + '_' + str(i+1) + '_on_' + str(i) + '.vcf'
                singleDonor = iterFold + "/donor_chr"+g+".fa"
                #singleHigh = iterFold + "/high_chr"+g+".vcf"
                singleHigh = iterFold + '/chr' + g + '_' + str(i+1) + '_high_confident_on_' + str(i) + '.vcf'
                #singleLow = iterFold + "/low_chr"+g+".vcf"
                singleLow = iterFold + '/chr' + g + '_' + str(i+1) + '_low_confident_on_' + str(i) + '.vcf'
                constructCommand = scriptPath + "/modify_genome_with_small_variants -r " + singleRef + " -s " + singleVcf + " -d " + singleDonor + " -h " + singleHigh + " -w " + singleLow
                single_command_list.append(constructCommand)
                recover_creating_command = args.python + ' ' + scriptPath + '/create_recover_position.py ' + singleHigh + '> ' + iterFold + '/chr' + g + '_' +  str(i+1) + '_to_' + str(i)
                single_command_list.append(recover_creating_command)
                #shell_run(recover_creating_command)
            
                if i > 0:
                    for target in range(i,0,-1):
                        previous_iter_fold = args.workspace + "/" + str(target)
                        recovering_all_command = args.python + ' ' + scriptPath + '/recover_coordinate_donor2ref.py -v ' + singleVcf + ' -r '+ previous_iter_fold + '/chr' + g + '_' + str(target) + '_to_' + str(target-1) + ' -f ' + iterFold + '/chr' + g + '_' + str(i+1) + '_on_' + str(target-1) + '.vcf'
                        recovering_confident_command = args.python + ' ' + scriptPath + '/recover_coordinate_donor2ref.py -v ' + singleHigh + ' -r '+ previous_iter_fold + '/chr' + g + '_' + str(target) + '_to_' + str(target-1) + ' -f ' + iterFold + '/chr' + g + '_' + str(i+1) + '_high_confident_on_' + str(target-1) + '.vcf'
                        single_command_list.append(recovering_all_command)
                        single_command_list.append(recovering_confident_command)
                pool.apply_async(serial_exe, (single_command_list, ))
                #shell_run(constructCommandt
            pool.close()
            pool.join()
            
            print ("finish constructing new chromosomes.\n")
            # combine donor genome
            combineCommand = "cat " + iterFold + "/donor_chr* > " + fastaFile
            shell_run(combineCommand)
            #delCommand = "rm " + iterFold + "/donor_chr*"
            #shell_run(delCommand)
            #delCommand = "rm " + iterFold + "/ref_chr*"
            #shell_run(delCommand)

if __name__ == '__main__':
    main()
