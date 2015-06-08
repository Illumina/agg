#!/usr/bin/env python 

import sys,subprocess,os,argparse,gzip

def get_sample_id(f_name):
    for line in gzip.open(f_name):
        line1 = line.split("\t")
        if line1[0] == "#CHROM":
            return line1[9].strip()


def check_binaries(args):
    for fname in [args.agg+"/bcftools-1.2/bcftools",args.agg+"/gvcftools-0.16/bin/extract_variants",args.agg+"/gvcftools-0.16/bin/get_called_regions",args.agg+"/bcftools-1.2/htslib-1.2.1/bgzip",args.agg+"/bcftools-1.2/htslib-1.2.1/tabix",args.agg+"/canon",args.agg+"/agg"]:
        if(not os.path.isfile(fname)):
            sys.exit(fname+" not found!\nCheck your provided paths.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='updates an agg database with new samples')
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
    for p in [db_dir+"/"+val for val in ["","variants","blocks","tmp","scripts"]]:
        if not os.path.isdir(p):
            sys.exit(p+" does not exist!")

    files = args.files
    if(not os.path.isfile(files)): sys.exit(files+"not found!")
    sample_list = {}
    
    sample_db = {}
    print "Checking original database..."
    ngood=0
    nsample=0
    for line in open("%s/samples.txt" % db_dir, "rb"):
        sampleid,v,b = line.split()
        variant_file=args.db+"/variants/"+sampleid+".bcf.csi"
        block_file=args.db+"/blocks/"+sampleid+".bed.gz.tbi"
        if not os.path.isfile(variant_file):
            sys.stderr.write("WARNING: %s does not exist. Will reprocess sample %s if gvcf was provided\n"%(variant_file,sampleid))
        elif not os.path.isfile(block_file):
            sys.stderr.write("WARNING: %s does not exist. Will reprocess sample %s if gvcf was provided\n"%(block_file,sampleid))
        else:                        
            sample_db[sampleid] = (v,b)
            ngood+=1
        nsample+=1
    print "There are",nsample-ngood,"/",nsample,"samples in the database that need re-ingesting..."

    sample_db_outfile = open("%s/samples.txt" % db_dir, "wb")
    for sampleid in sample_db:
        variant_file=args.db+"/variants/"+sampleid+".bcf"
        bed_file=args.db+"/blocks/"+sampleid+".bed.gz"
        sample_db_outfile.write("\t".join([sampleid, variant_file, bed_file]))
        sample_db_outfile.write("\n")

    jids = []
    ingestion_list = [(get_sample_id(val.strip()),val.strip()) for val in open(files, "rb").read().split()]
    ingestion_list = [(sampleid,gvcf) for sampleid,gvcf in ingestion_list if sampleid not in sample_db]
    for sampleid,gvcf in ingestion_list:     
        print sampleid,gvcf

    for i in xrange(0,len(ingestion_list),args.nprocess):
        chunk = [val for val in ingestion_list[i:(i+args.nprocess)]]
        input_gvcfs = [gvcf for sampleid,gvcf in chunk]
        for sampleid,gvcf in chunk:
            if not os.path.isfile(gvcf):
                sys.stderr.write("Warning: %s does not exist!\n"%gvcf)
            else:
                variant_file="%s/variants/%s.bcf" %(db_dir, sampleid)
                bed_file="%s/blocks/%s.bed.gz" % (db_dir, sampleid)
                sample_db_outfile.write("\t".join([sampleid, variant_file, bed_file]))
                sample_db_outfile.write("\n")

        if args.batch=="qsub":
            qsub = ["qsub","-pe","threaded",str(args.nprocess),"-e",db_dir+"/tmp","-o",db_dir+"/tmp","-terse","-S","/bin/bash","-cwd","-b","y",sys.executable,scripts+"ingest_multi_gvcf.py"]
        else: 
            qsub = ["sbatch","-V","--time=1:00:00","--parsable",scripts+"ingest_multi_gvcf.py"]

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
        cmd = ["sbatch","-V","--time=2:00:00","-d",hold,"--parsable",scripts+"multi_collate.py",db_dir,"-agg",args.agg+"/","-nprocess",str(args.nprocess)]
    print " ".join(cmd)

    subprocess.call(cmd)
    print "All jobs submitted.\nPlease wait until jobs finish before querying database %s"%db_dir


