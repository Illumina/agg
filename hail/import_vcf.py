import argparse,time,sys,os
import pandas as pd, numpy as np

MAXDEPTH = 4

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='hail script to apply a simple set of hard filters')
    parser.add_argument('-files', metavar='files', type=str, help='file containing list of vcfs to import',required=True)
    parser.add_argument('-vds', metavar='vds', type=str, help='vds prefix',required=True)
    parser.add_argument('-lcr', metavar='lcr', type=str, help='bed file containing the low-complexity regions to filter',default=None)
    parser.add_argument('-tmp', metavar='tmp', type=str, help='tmp directory',default="/tmp/")
    parser.add_argument('-rename', metavar='rename', type=str, help='two column text file for sample renaming',default=None)
    parser.add_argument('--vcf',action='store_true',default=False,help='output a big vcf')
     
    args = parser.parse_args()

    import hail,os,collections

    print "VCF:",args.files
    print "VDS:",args.vds
    print "LCR:",args.lcr
    
    if os.path.exists(args.vds):
        print "ERROR: ",args.vds,"exists!"
        sys.exit()
            
    hc = hail.HailContext(log="hail.log",tmp_dir=args.tmp)

    time0 = time.time()
    vds=hc.import_vcf(open(args.files).read().strip().split(),force_bgz=True,store_gq=True).split_multi().variant_qc()
    vds=vds.annotate_variants_expr(['va.filters = [""][:0].toSet()',"va.info.AC_raw = va.qc.AC"]).annotate_variants_expr('va.ft.alt_gq_median = gs.filter(g=>g.isCalledNonRef).map(g=>g.gq).collect().median()')
    print "VCF conversion took",time.time()-time0,"seconds"

    if args.lcr!=None:
        time0 = time.time()
        vds = vds.annotate_variants_table(hail.KeyTable.import_bed(args.lcr), root='va.lowComplexityRegion').annotate_variants_expr('va.filters = if(va.lowComplexityRegion) va.filters.add("LCR") else va.filters')
        print "LCR annotation took",time.time()-time0,"seconds"                                      

    if args.rename!=None:
        sample_rename = dict([val.strip().split() for val in open(args.rename)])
        vds = vds.rename_samples(sample_rename)

    time0 = time.time()        
    vds=vds.filter_genotypes("g.gq>=20 && g.dp>=10 && (!g.isHet() || (g.ad[1]/g.ad.sum())>=0.2) ").variant_qc().annotate_variants_expr(["va.info.InbreedingCoeff = if(va.qc.AC>0) (1-va.qc.nHet/(va.qc.nCalled*2*va.qc.AF*(1-va.qc.AF))) else 0.0","va.info.AC = va.qc.AC","va.info.AF = va.qc.AF"])

    filter_expressions = ['va.filters = if(!isMissing(va.info.InbreedingCoeff) && va.info.InbreedingCoeff < -0.3) va.filters.add("InbreedingCoeff") else va.filters',
                          'va.filters = if(va.info.AC == 0) va.filters.add("AC0") else va.filters',
                          'va.filters = if(!isMissing(va.ft.alt_gq_median) && va.ft.alt_gq_median<20) va.filters.add("LOWGQ") else va.filters',
                          'va.filters = if(va.qc.nCalled<0.9) va.filters.add("LOWCALL") else va.filters']

    for e in filter_expressions:
        vds=vds.annotate_variants_expr(e)
#    vds=vds.annotate_variants_expr(filter_expressions) ##this does not work (i think because the expressions are modifying the same field)
    vds=vds.annotate_variants_expr('va.pass = va.filters.isEmpty()')
    print "Filter pass one took ",time.time()-time0,"seconds"
    
    time0 = time.time()        
    ##sets up a simple max depth filter
    sample_expressions=['sa.altDepthStats = gs.filter(g => va.info.AC>=10 && va.pass && g.isCalledNonRef() && g.dp<1000).map(g=>g.dp).stats()']
    vds=vds.annotate_samples_expr(sample_expressions)
    vds=vds.annotate_variants_expr('va.ft.alt_dp_mean = gs.filter(g=>g.isCalledNonRef).map(g=>(g.dp-sa.altDepthStats.mean)/sa.altDepthStats.stdev).stats().mean')
    vds=vds.annotate_variants_expr('va.filters = if(!isMissing(va.ft.alt_dp_mean) && va.ft.alt_dp_mean>%f) va.filters.add("HIGHDP") else va.filters'%MAXDEPTH)
    vds=vds.annotate_variants_expr('va.pass = va.filters.isEmpty()')
    print "Filter pass two took ",time.time()-time0,"seconds"

    vds = vds.set_va_attributes('va.filters', {'AC0': 'no alternate genotypes passed per-genotype hard filters',
                                               'LCR': 'variant falls in a low-complexity region',
                                               'InbreedingCoeff': 'inbreeding coefficient < -0.3 (excessive heterozygosity)',
                                               'HIGHDP': 'alternate genotypes have excessively high depth',
                                               'LOWGQ': 'the median GQ at alternate genotypes was <20',
                                               'LOWCALL':'<0.9 genotypes had a high quality genotype call'})    

    time0=   time.time()
    vds.write(args.vds)   
    print "VDS write took",time.time()-time0,"seconds"                                      
    variantqc_table = vds.variants_table().to_pandas()

    if args.vcf:
        vds.export_vcf(args.vds+"vcf.bgz",parallel=True)
    raw_counts= vds.count()

    vds_pass = vds.filter_variants_expr('va.pass').sample_qc()
    pass_counts = vds_pass.count()    
    ## some simply summaries of variants (I think this is inefficient)
    s = vds_pass.samples_table().to_pandas()
    s.to_csv("samples.csv")
    print "\n\nSample QC (PASS variants):\n"
    tab = pd.DataFrame({'nSNP':s['sa.qc.nSNP'].describe(),'TiTv':s['sa.qc.rTiTv'].describe(),'nIndel':(s['sa.qc.nInsertion']+s['sa.qc.nDeletion']).describe(),"nSingleton":s['sa.qc.nSingleton'].describe()}).round(2)
    print tab.ix[[1,3,4,5,6,7]]

    print "\n\nCohort wide summary:"
    q1=['variants.filter(v => va.info.AC==1 && v.altAllele().isSNP()).count()',
        'variants.filter(v => va.info.AC==1 && v.altAllele().isIndel()).count()',
        'variants.filter(v => va.info.AC==1 && v.altAllele().isSNP() && va.pass).count()',
        'variants.filter(v => va.info.AC==1 && v.altAllele().isIndel() && va.pass).count()']
    q2=['variants.filter(v => va.info.AF<0.01 && va.info.AC>1 && v.altAllele().isSNP()).count()',
        'variants.filter(v => va.info.AF<0.01 && va.info.AC>1 && v.altAllele().isIndel()).count()',
        'variants.filter(v => va.info.AF<0.01 && va.info.AC>1 && v.altAllele().isSNP() && va.pass).count()',
        'variants.filter(v => va.info.AF<0.01 && va.info.AC>1 && v.altAllele().isIndel() && va.pass).count()']
    q3=['variants.filter(v => va.info.AF>=0.01 && v.altAllele().isSNP()).count()',
        'variants.filter(v => va.info.AF>=0.01 && v.altAllele().isIndel()).count()',
        'variants.filter(v => va.info.AF>=0.01 && v.altAllele().isSNP() && va.pass).count()',
        'variants.filter(v => va.info.AF>=0.01 && v.altAllele().isIndel() && va.pass).count()']
    q4=['variants.filter(v => v.altAllele().isSNP()).count()',
        'variants.filter(v => v.altAllele().isIndel()).count()',
        'variants.filter(v => v.altAllele().isSNP() && va.pass).count()',
        'variants.filter(v => v.altAllele().isIndel() && va.pass).count()'   ]

    
    counts = vds.query_variants(q1 + q2 + q3 + q4)

    print pd.DataFrame([['singleton']+counts[:4],['rare (AF<.01)']+counts[4:8],['common (AF>=.01)']+counts[8:12],['all']+counts[12:16]],columns=["frequency","snp.all","indel.all","snp.pass","indel.pass"])
    

