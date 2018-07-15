import argparse
import gzip
import os
import subprocess
import shutil, sys
from subprocess import Popen, PIPE

def main():
    args = parse_args()
    bedPath=args.bedFolder
    bamAPath=args.bamFolderA
    bamBPath=args.bamFolderB
    resultPath=args.resultFolder
    FolderName = args.AorB

    print(FolderName)

    ## select source and Destination files based on A available in B or opposit
    if FolderName == 'A':
        bamPathsource = bamAPath
        bamPathDestenation = bamBPath
        resultfolder = resultPath + "/A"
        bedFilePath = bedPath+"/A_SNP.bed"
    else:
        bamPathsource = bamBPath
        bamPathDestenation = bamAPath
        resultfolder = resultPath + "/B"
        bedFilePath = bedPath+"/B_SNP.bed"

    #shutil.rmtree(resultPath)
    #if not os.path.exists(resultPath):
       # os.makedirs(resultPath)
    #os.makedirs(resultfolder)

    ## create Availabalities files
    get_GenomePositionAvailability_fromX_inX(bamPathsource, bamPathDestenation, resultfolder, bedFilePath)
    print("Done.")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-bed", "--bedFolder")
    parser.add_argument("-bamA", "--bamFolderA")
    parser.add_argument("-bamB", "--bamFolderB")
    parser.add_argument("-result", "--resultFolder")
    parser.add_argument("-AorB", "--AorB")
    args = parser.parse_args()
    return parser.parse_args()

def get_genome_position_range(line):
    line = line.split("\n")[0]
    return line.split("\t")[0]+":"+line.split("\t")[1]+"-"+line.split("\t")[1]


def Get_Set(bamFile):
    hashset = set()
    with open(bamFile,"rb") as f:
        for line in f:
            hashset.add(line.rstrip('\t\n\r'))
    return hashset

def get_ID(bamFile, Genome_position):
    genomeIDs = "samtools view "+bamFile+" "+Genome_position+" | cut -f 1 |rev |cut -f -4 -d ':' |rev| cut -f 1 -d '#'"
    process = Popen(genomeIDs, shell=True, stdout=PIPE,stderr=PIPE)
    stdout, stderr = process.communicate()
    hashset = set()
    if len(stdout)>1:
        for ID in stdout.split("\n"):
            hashset.add(ID)
    return hashset



## for each position of VCF file for ID related to that position fron A bam and calculate the SUM/Count and flage of the appiarance of all ID's in B bam file
## the resulr will be in the csvPath folder for each sample
def get_GenomePositionAvailability_fromX_inX(bamPathsource,bamPathDestenation,resultfolder,bedFilePath):
    for bamFilename in os.listdir(bamPathsource):
        if bamFilename.endswith(".bam"):
            bam_Set = Get_Set(bamPathDestenation+"/"+bamFilename+".txt")
            with open(bedFilePath) as bedFile:
                for bed_line in bedFile:
                    if "#" not in bed_line:
                        Genome_position = get_genome_position_range(bed_line)
                        ID_set = get_ID(bamPathsource+"/"+bamFilename, Genome_position)
                        ID_Count = len(ID_set)
                        Match_Count = 0
                        if ID_Count > 0:
                            ID_Intersect = bam_Set & ID_set
                            Match_Count = len(ID_Intersect)
                            fraction= 0 if ID_Count == 0 else (Match_Count / float(ID_Count))
                        with open(resultfolder + "/" + bamFilename, 'a') as resultFilename:
                            resultFilename.write(Genome_position.split(":")[0]+"\t"+Genome_position.split("-")[1]+"\t"+str(fraction)+"\t"+str(Match_Count)+"\t"+str(ID_Count)+"\n")


if __name__ == '__main__': main()