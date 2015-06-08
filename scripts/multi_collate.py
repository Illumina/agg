#!/usr/bin/env python 

import sys,subprocess,os,argparse,time,gzip,multiprocessing,re

def readChunks(bcftools,fname):
    chunk_size = 30000000
    contig_info=subprocess.check_output("%s view -h %s | grep contig"%(bcftools,fname),shell=True).split("\n")
    reg = []
    for line in contig_info:
        if "##contig=<ID=" in line:
            chrom = (line.split("ID=")[1].split(",")[0])
            if "M" not in chrom:
                start = 0
                stop = int(line.split("length=")[1].split(">")[0])
                reg+=["%s:%d-%d"%(chrom,i,i+chunk_size-1) for i in range(int(start),int(stop),chunk_size)]

    return(reg)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='ingests a gvcf into an agg database')
    parser.add_argument('db', metavar='db', type=str, help='database dir')    
    parser.add_argument('-agg', metavar='agg', type=str, help='path to agg install',required=True)
    parser.add_argument('-nprocess', metavar='nprocess', type=int, help='number of processes to use',default=16)
    args = parser.parse_args()

    if not os.path.isdir(args.db+"/variants/"):
        sys.exit(args.db+"/variants does not exist!")

    if not os.path.isdir(args.db+"/blocks/"):
        sys.exit(args.db+"/blocks does not exist!")

    ##binaries
    bcftools = args.agg+"/bcftools-1.2/bcftools"
    agg = args.agg+"/agg"


    vcf = args.db+"/variants/"+open(args.db+"/samples.txt").next().split()[0] + ".bcf"
    chunks = readChunks(bcftools,vcf)
    if(len(chunks)==0):
        sys.exit("ERROR: problem creating regions to process in file %s"%vcf)

    ##collate chunks with multiprocessing pool
    def collate_chunk(region):
        tmp1 = "%s/tmp/%s.1.bcf"%(args.db,region)
        tmp2 = "%s/tmp/%s.bcf"%(args.db,region)
        subprocess.call("%s/agg collate %s %s %s > %s.stdout 2> %s.stderr"%(args.agg,args.db,tmp1,region,tmp1,tmp1),shell=True)
        if not os.path.isfile(tmp1):
            raise Exception("error running collate on region "+region+". Problem with ingestion? Check "+tmp1+".stderr for details")
        subprocess.call("%s/canon %s > %s "%(args.agg,tmp1,tmp2),shell=True)
        if not os.path.isfile(tmp2):
            raise Exception("error running canon on region "+region+". Problem with collation?")            
        return(tmp2)

    print "Collating variants..."
    pool = multiprocessing.Pool(processes=args.nprocess)
    output_files = pool.map(collate_chunk, chunks)

    for f in output_files:
        if not os.path.isfile(f):
            sys.exit("Problem creating %s"%f)

    ##concatenate the final chunks into our database
    print "Concatenating..."
    cmd = bcftools + " concat " + " ".join(output_files) + " -O b -o "+args.db+"/sites.bcf"
    print cmd
    subprocess.call(cmd,shell=True)
    subprocess.call(bcftools + " index "+args.db+"/sites.bcf",shell=True)
    print "Cleaning up..."
    for f in output_files:
        os.remove(f)

    
