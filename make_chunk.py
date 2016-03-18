#!/usr/bin/env python 

import sys,subprocess,os,argparse,time,multiprocessing,tempfile,shutil


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='makes an agg chunk')
    parser.add_argument('input', metavar='input', type=str, help='plain text file listing gvcfs to be ingested')
    parser.add_argument('-output', metavar='output', type=str, help='output prefix',required=True)
    parser.add_argument('-ref', metavar='ref', type=str, help='reference genome',required=True)
    parser.add_argument('-tmp', metavar='tmp', type=str, help='tmp directory',default="/tmp/")
    parser.add_argument('-nprocess', metavar='nprocess', type=int, help='number of processes to use (defaults to number of CPUs found)',default=None)
    parser.add_argument('-agg', metavar='agg', type=str,default=None, help='agg binary (defaults to binary in script dir)')
    parser.add_argument('--allow-duplicates', dest='allow_duplicates', action='store_const', const=True,default=False, help='allow duplicate sample names (not recommended)')
    parser.add_argument('--ignore-non-matching-ref', action='store_true',default=False)
    args = parser.parse_args()
    sys.stderr.write(" ".join(sys.argv)+"\n")

    ##binaries
    if args.agg==None:
        args.agg=os.path.dirname(os.path.realpath(sys.argv[0]))+"/agg"

    if not os.path.isfile(args.agg):
        sys.exit(args.agg + " binary not found!")

    if os.path.isfile(args.output):
        sys.exit(args.output+" already exists! will not overwrite")
    
    if not os.path.isdir(os.path.dirname(args.output)):
        sys.exit("Output directory "+os.path.abspath(args.output)+" does not exist")

    if not os.path.isfile(args.input):
        sys.exit(args.input+" does not exist!")

    if args.nprocess==None:
        args.nprocess=multiprocessing.cpu_count()
    sys.stderr.write("Using %d processes\n"%args.nprocess)
    args.tmp+="/"
    tmp_dir = tempfile.mkdtemp(prefix=args.tmp)
    sys.stderr.write("tmp dir: %s\n"%tmp_dir)
    assert os.path.isdir(tmp_dir)
    gvcfs = open(args.input).read().strip().split()
    sys.stderr.write("processing %d gvcfs\n"%len(gvcfs))

    ##collate chunks with multiprocessing pool
    def process_gvcf(f):  
        tmp_out=tempfile.mktemp(prefix=os.path.basename(f).split(".")[0],dir=tmp_dir)
#        tmp_out = "%s/%s"%(tmp_dir,os.path.basename(f).split(".")[0])
        cmd = args.agg + " ingest1 " + f + " -o " + tmp_out + " -f " + args.ref
        if args.ignore_non_matching_ref:
            cmd += " --ignore-non-matching-ref"
        sys.stderr.write(cmd+"\n")
        print cmd
        try:
            subprocess.check_output(cmd,shell=True)
            return(tmp_out+".bcf")
        except:
            sys.stderr.write("WARNING: Problem processing %s. Will not be included"%f)
            return None

    try:
        time0=time.time()
        sys.stderr.write("Running agg ingest1...\n")
        pool = multiprocessing.Pool(processes=args.nprocess)
        output_files = pool.map(process_gvcf, gvcfs)
        sys.stderr.write("ingest1 took %f seconds\n"%(time.time()-time0))
        output_files.sort()
        fout = open("%s/ingest1.txt"%tmp_dir,"wb")
        fout.write("\n".join([val for val in output_files if val!=None]))
        del fout

        ##concatenate the final chunks into our database
        time0=time.time()
        sys.stderr.write("running agg ingest2...\n")
        cmd = "%s ingest2 -l %s -@%d -o %s"%(args.agg,"%s/ingest1.txt"%tmp_dir,args.nprocess,args.output)
        if args.allow_duplicates:
            cmd += " --allow-duplicates"
        sys.stderr.write(cmd+"\n")
        subprocess.check_output(cmd,shell=True)
        sys.stderr.write("ingest2 took %f seconds\n"%(time.time()-time0))
#        shutil.rmtree(tmp_dir)
    except:
        sys.stderr.write("there was a problem. removing temporary files and exiting\n\n")
#        shutil.rmtree(tmp_dir)
        sys.exit(1)
