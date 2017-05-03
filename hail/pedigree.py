import argparse,time,sys,tempfile
import hail,os,collections
import pandas as pd, numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='hail script to perform a few basic pedigree analyses')
    parser.add_argument('-vds', metavar='vds', type=str, help='vds prefix',required=True)
    parser.add_argument('-tmp', metavar='tmp', type=str, help='tmp directory',default="/tmp/")
    parser.add_argument('-pedigree', metavar='pedigree', type=str, help='plink pedigree file (will remove samples not in the file)',required=True)
    parser.add_argument('-gnomad', metavar='gnomad', type=str, help='gnomad .vds to annotate allele frequencies from',required=True)    

    args = parser.parse_args()



    print "VDS:",args.vds
    print "Pedigree:",args.pedigree
                        
    hc = hail.HailContext(log="/dev/null")

    time0 = time.time()
    vds = hc.read(args.vds)

    ids_to_keep = [val.split()[1] for val in open(args.pedigree)]
    vds=vds.filter_samples_list(ids_to_keep).annotate_samples_fam(args.pedigree)


    tmp_dir = tempfile.mkdtemp(prefix=args.tmp+"/")
    vds.mendel_errors(tmp_dir+"/1", args.pedigree)
    mendel = pd.read_csv(tmp_dir+"/1.fmendel",delim_whitespace=True)         

    print "Problematic pedigrees:"
    print  mendel[mendel.NSNP > (mendel['NSNP'].median() + 6* mendel['NSNP'].mad())]
    print "Parental relationships will be deleted"
        
    invalid_parents = mendel.MAT[mendel.NSNP > (mendel['NSNP'].median() + 6* mendel['NSNP'].mad())].tolist() + mendel.PAT[mendel.NSNP > (mendel['NSNP'].median() + 6* mendel['NSNP'].mad())].tolist()
    
    fam = pd.read_csv(args.pedigree,delim_whitespace=True,names=["FID","IID","DAD","MUM","SEX","PHENO"])
    fam = fam.set_value(np.logical_or(np.in1d(fam.DAD,invalid_parents),np.in1d(fam.MUM,invalid_parents)),"DAD","0")
    fam = fam.set_value(np.logical_or(np.in1d(fam.DAD,invalid_parents),np.in1d(fam.MUM,invalid_parents)),"MUM","0")

    fam.to_csv(tmp_dir+"/corrected.fam",sep=" ",header=False,index=False)
    newped=tmp_dir+"/corrected.fam"
#    print vds.variant_schema
    vds=vds.annotate_samples_fam(newped).tdt(newped).annotate_samples_expr('sa.founder = isMissing(sa.fam.patID) && isMissing(sa.fam.patID)').variant_qc().filter_variants_expr("va.pass && va.qc.AC>0").annotate_variants_expr('va.founder = gs.filter(g => sa.founder).callStats(g => v)')
                        
    variantqc_table = vds.variants_keytable().to_pandas()
    variantqc_table.columns = [val.replace(".","_") for val in variantqc_table.columns]

    variantqc_table['founder_ac'] = [val[1] for val in variantqc_table['va_founder_AC']]
    variantqc_table['founder_af'] = variantqc_table.founder_ac / variantqc_table.va_founder_AN
    variantqc_table['is_snp'] = [(len(val[0][0])==1 and len(val[0][1])==1) for val in variantqc_table['v_altAlleles']]
    variantqc_table['is_snp'] = [(len(val[0][0])==1 and len(val[0][1])==1) for val in variantqc_table['v_altAlleles']]                        
#    variantqc_table.to_csv("variantqc_summary.csv")
    
    vtmp = variantqc_table.query("is_snp & founder_ac==1")
    print "Singleton SNP transmission rate       = %f"%(float(vtmp.va_tdt_nTransmitted.sum())/(vtmp.va_tdt_nTransmitted+vtmp.va_tdt_nUntransmitted).sum())
    vtmp = variantqc_table.query("not is_snp & founder_ac==1")
    print "Singleton indel transmission rate     = %f"%(float(vtmp.va_tdt_nTransmitted.sum())/(vtmp.va_tdt_nTransmitted+vtmp.va_tdt_nUntransmitted).sum())
    vtmp=variantqc_table.query("is_snp & founder_ac>1")
    print "Non-singleton SNP transmission rate   = %f"%(float(vtmp.va_tdt_nTransmitted.sum())/(vtmp.va_tdt_nTransmitted+vtmp.va_tdt_nUntransmitted).sum())
    vtmp=variantqc_table.query("not is_snp & founder_ac>1")
    print "Non-singleton indel transmission rate = %f"%(float(vtmp.va_tdt_nTransmitted.sum())/(vtmp.va_tdt_nTransmitted+vtmp.va_tdt_nUntransmitted).sum())
