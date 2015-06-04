#!/usr/bin/env python 

import sys,subprocess,os,argparse,time,gzip,multiprocessing

def get_sample_id(f_name):
    for line in gzip.open(f_name):
        line1 = line.split("\t")
        if line1[0] == "#CHROM":
            return line1[9].strip()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='ingests a gvcf into an agg database')
    parser.add_argument('gvcf', metavar='gvcf', type=str, help='gvcf')
    parser.add_argument('db', metavar='db', type=str, help='database dir')    
    parser.add_argument('-agg', metavar='agg', type=str, help='path to agg install',required=True)
    parser.add_argument('-ref', metavar='ref', type=str, help='path to reference genome',required=True)
    args = parser.parse_args()

    if not os.path.isdir(args.db+"/variants/"):
        sys.exit(args.db+" does not exist!")

    if not os.path.isdir(args.db+"/blocks/"):
        sys.exit(args.db+" does not exist!")

    time0 = time.time()

    os.environ["BCFTOOLS_PLUGINS"] = args.agg+"/bcftools-1.2/plugins/"
    os.environ["LD_LIBRARY_PATH"] = args.agg+"/bcftools-1.2/htslib-1.2.1/"
    sampleid=get_sample_id(args.gvcf)

    ##binaries
    bcftools = args.agg+"/bcftools-1.2/bcftools"
    get_called_regions = args.agg+"/gvcftools-0.16/bin/get_called_regions"
    bgzip = args.agg+"/bcftools-1.2/htslib-1.2.1/bgzip"
    extract_variants = args.agg+"/gvcftools-0.16/bin/extract_variants"
    canon = args.agg+"/canon"
    tabix = args.agg+"/bcftools-1.2/htslib-1.2.1/tabix"

    ##output names
    outbcf = "%s/variants/%s.bcf"%(args.db,sampleid)
    outbed = "%s/blocks/%s.bed.gz"%(args.db,sampleid)
    print("Outputting:")
    print(outbcf)
    print(outbed)

    if  os.path.isfile(outbcf):
        sys.exit("%s already exists!"%outbcf)

    if  os.path.isfile(outbed):
        sys.exit("%s already exists!"%outbed)
    
    print("Making %s..."%outbcf)
#    norm = [bcftools,"norm","-m","-any","|",bcftools,"+fixploidy","|",bcftools,"norm","-m","+any","-","-O","u","|",bcftools,"norm","-m","-any","-","-O","u"]
    norm = [bcftools,"norm","-m","-any","|",bcftools,"norm","-f",args.ref,"-","-O","u","|",bcftools,"+fixploidy","--","-p","/dev/null","|",bcftools,"norm","-m","+any","-","-O","u","|",bcftools,"norm","-m","-any","-","-O","u"] 
    cmd=["zcat",args.gvcf,"|",extract_variants,"|",bcftools,"annotate","-x","INFO","-O","u","|"] + norm + ["|",canon,"-",">",outbcf]
    print( " ".join(cmd) )
    subprocess.call(" ".join(cmd),shell=True)
    subprocess.call([bcftools,"index",outbcf])

    print("Making %s..."%outbed)
    cmd = ["zcat",args.gvcf,"|",get_called_regions,"|",bgzip,"-c",">",outbed]
    print( " ".join(cmd))
    subprocess.call(" ".join(cmd),shell=True)
    subprocess.call([tabix,outbed])

    print("Ingestion of %s took %f seconds"%(args.gvcf,time.time()-time0))
