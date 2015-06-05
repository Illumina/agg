#!/usr/bin/env python 

import sys,subprocess,os,argparse,time,gzip,multiprocessing

def get_sample_id(f_name):
    for line in gzip.open(f_name):
        line1 = line.split("\t")
        if line1[0] == "#CHROM":
            return line1[9].strip()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='ingests a gvcf into an agg database')
    parser.add_argument('gvcf', metavar='gvcf', type=str, help='comma separated list gvcfs (1 process is started per gvcf)')
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


    ##binaries
    vt  = args.agg+"/vt-0.57/vt"
    bcftools = args.agg+"/bcftools-1.2/bcftools"
    get_called_regions = args.agg+"/gvcftools-0.16/bin/get_called_regions"
    bgzip = args.agg+"/htslib-1.2.1/bgzip"
    extract_variants = args.agg+"/gvcftools-0.16/bin/extract_variants"
    canon = args.agg+"/canon"
    tabix = args.agg+"/htslib-1.2.1/tabix"

    ##sub command that will be called per gvcf
    def ingest(fname):
        ##output names
        sampleid=get_sample_id(fname)
        outbcf = "%s/variants/%s.bcf"%(args.db,sampleid)
        outbed = "%s/blocks/%s.bed.gz"%(args.db,sampleid)
        print("Outputting:")
        print(outbcf)
        print(outbed)

        # if  os.path.isfile(outbcf):            sys.exit("%s already exists!"%outbcf)
        # if  os.path.isfile(outbed):            sys.exit("%s already exists!"%outbed)

        print("Making %s..."%outbcf)
        subprocess.call('%s view -h %s | %s annotate  -x ^INFO/SNVSB  |  sed s/"##FORMAT=<ID=AD,Number=."/"##FORMAT=<ID=AD,Number=R"/g > %s.hdr'%(bcftools,fname,bcftools,outbcf),shell=True)

        norm = "%s annotate -Ou -x ^INFO/SNVSB | %s reheader -h %s.hdr - | %s +fixploidy -Ou | %s norm -m -any | %s normalize -r %s - "%(bcftools,bcftools,outbcf,bcftools,bcftools,vt,args.ref)

        cmd=["zcat",fname,"|",extract_variants,"|"] + [norm] + ["|",canon,"-",">",outbcf]

        print( " ".join(cmd) )

        subprocess.call(" ".join(cmd),shell=True)
        subprocess.call([bcftools,"index",outbcf])
        os.remove("%s.hdr"%outbcf)
        print("Making %s..."%outbed)
        cmd = ["zcat",fname,"|",get_called_regions,"|",bgzip,"-c",">",outbed]
        print( " ".join(cmd))
        subprocess.call(" ".join(cmd),shell=True)
        subprocess.call([tabix,outbed])

        print("Ingestion of %s took %f seconds"%(fname,time.time()-time0))

    gvcf_list = args.gvcf.split(",")
    print gvcf_list
    nfile = len(gvcf_list)
    pool = multiprocessing.Pool(processes=nfile)
    pool.map(ingest, gvcf_list)

