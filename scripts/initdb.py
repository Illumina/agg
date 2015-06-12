#!/usr/bin/env python 

import sys,subprocess,os,argparse,gzip

def get_sample_id(f_name):
    for line in gzip.open(f_name):
        line1 = line.split("\t")
        if line1[0] == "#CHROM":
            return line1[9].strip()


def check_binaries(args):
    vt  = args.agg+"/vt-0.57/vt"
    bcftools = args.agg+"/bcftools-1.2/bcftools"
    get_called_regions = args.agg+"/gvcftools-0.16/bin/get_called_regions"
    bgzip = args.agg+"/htslib-1.2.1/bgzip"
    extract_variants = args.agg+"/gvcftools-0.16/bin/extract_variants"
    canon = args.agg+"/canon"
    tabix = args.agg+"/htslib-1.2.1/tabix"
    for f in [vt,bcftools,get_called_regions,bgzip,extract_variants,canon,tabix]:
        if not os.path.isfile(f):
            sys.exit(f+" not found!\nCheck your provided paths.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='initialises an agg database')
    parser.add_argument('-db', metavar='db', type=str, help='database to be created',required=True)
    parser.add_argument('-files', metavar='files', type=str, help='plain text file listing gvcfs to be ingested',required=True)
    parser.add_argument('-agg', metavar='agg', type=str, help='path to agg install',required=True)
    parser.add_argument('-chrom_prefix', metavar='chrom_prefix', type=str, help='prefix for chromosomes',default="")
    parser.add_argument('-ref', metavar='ref', type=str, help='path to reference genome',required=True)
    parser.add_argument('-batch', metavar='batch', type=str, help='queue submission command [qsub,sbatch]',default="qsub")
    parser.add_argument('-nprocess', metavar='nprocess', type=int, help='number of ingestion jobs per node',default=1)
    args = parser.parse_args()

    check_binaries(args)
    if args.batch not in ["qsub","sbatch"]:
        sys.exit("ERROR: unsupported queueing system: " + args.batch)

    scripts =  os.path.dirname(os.path.realpath(__file__)) +"/"
    paths = ["-agg",args.agg,"-ref",args.ref]
    bcftools = args.agg+"/bcftools-1.2/bcftools"


    db_dir = args.db
    files = args.files
    if(not os.path.isfile(files)): sys.exit(files+"not found!")
    sample_list = {}
    if not os.path.isdir(db_dir):
        subprocess.call(["mkdir", "-p", db_dir + "/variants"])
        subprocess.call(["mkdir", "-p", db_dir + "/blocks"])
        subprocess.call(["mkdir", "-p", db_dir + "/scripts"])
        subprocess.call(["mkdir", "-p", db_dir + "/tmp"])
    else:
        sys.exit(db_dir + " already exists!")

    sample_db = open("%s/samples.txt" % db_dir, "wb")
    jids = []
    gvcfs = [val.strip() for val in open(files, "rb").read().split()]

    for i in xrange(0,len(gvcfs),args.nprocess):
        input_gvcfs= gvcfs[i:(i+args.nprocess)]

        for gvcf in input_gvcfs:
            if not os.path.isfile(gvcf):
                sys.stderr.write("Warning: %s does not exist!\n"%gvcf)
            else:
                sample_id = get_sample_id(gvcf)
                variant_file="%s/variants/%s.bcf" %(db_dir, sample_id)
                bed_file="%s/blocks/%s.bed.gz" % (db_dir, sample_id)
                sample_db.write("\t".join([sample_id, variant_file, bed_file]))
                sample_db.write("\n")


        if args.batch=="qsub":
            qsub = ["qsub","-q","devel-s.q","-pe","threaded",str(args.nprocess),"-e",db_dir+"/tmp","-o",db_dir+"/tmp","-terse","-S","/bin/bash","-cwd","-b","y",sys.executable,scripts+"ingest_multi_gvcf.py"]
        else: 
            qsub = ["sbatch","--reservation=GEL","--time=1:00:00","--parsable",scripts+"ingest_multi_gvcf.py"]

        cmd = qsub + [",".join(input_gvcfs), db_dir] + paths        
        print " ".join(cmd)
        jids.append(subprocess.check_output(cmd).strip())

    if len(jids)==0:
        sys.exit("ERROR: no jobs submitted")

    hold = ",".join(jids)
    jids = []

    if args.batch=="qsub":        
        cmd = ["qsub","-pe","threaded",str(args.nprocess),"-e",db_dir+"/tmp","-o",db_dir+"/tmp","-hold_jid",hold,"-terse","-N","collate","-S","/bin/bash","-cwd","-b","y",sys.executable,scripts+"multi_collate.py",db_dir,"-agg",args.agg+"/","-nprocess",str(args.nprocess)]
    else: 
        cmd = ["sbatch","--reservation=GEL","--time=2:00:00","-d",hold,"--parsable",scripts+"multi_collate.py",db_dir,"-agg",args.agg+"/","-nprocess",str(args.nprocess)]
    print " ".join(cmd)

    subprocess.call(cmd)
    print "All jobs submitted.\nPlease wait until jobs finish before querying database %s"%db_dir


