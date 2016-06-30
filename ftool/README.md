# filter_tools
Population based variant filtering pipeline. Requires a (preferably site-only) vcf from agg.

Build:
```
git clone git@git.illumina.com:joconnell/filter_tools.git
cd filter_tools/
make
./ftool
```
Usage:
```
$ ./ftool 

Program:	ftool (some post-hoc filtering for agg output)
Contact:	joconnell@illumina.com

Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.

Usage:	ftool <command> [options]

Commands:

annotate        adds some INFO fields that are useful for filtering
evaluate        tabulates mendel rates/tstv across quantiles of an INFO field
power           looks at proportion of "true" variants maintained under filtering routines
```

##ftool annotate
Calculates and annotates useful annotations for variant filtering. The most useful thing it adds at the moment is `INFO/S_DP` which is an "frequency corrected depth measure", `S_DP` is the number of standard deviations away from the expectation. [See these slides for a summary](https://confluence.illumina.com/download/attachments/108066684/20160308.agg.filtering.pptx?version=1&modificationDate=1459255893050&api=v2).

```
ftool annotate -i 'QUAL>=20 && TYPE=="snp" && DP<200000 && DPF<DPF' sites.vcf.gz | bgzip -c > sites.annotate.vcf.gz
```
The model does two passes
1. calculate the median/MAD for depth in different allele frequency bins on the *filtered* sites
2. calcualte the standard deviation for each DP value based on the estimate parameters on *all* sites

The reason for the `-i` filtering is we want to calculate parameters on some reasonably sensible DP values, even though we are using robust estimators, some light filtering is desirable.

