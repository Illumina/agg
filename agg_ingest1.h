#include "agg.h"
#include "kseq.h"
#include "kstring.h"

KSTREAM_INIT(gzFile, gzread, 16384)
#include "vcfnorm.h"

extern "C" {
#include "htslib/hts.h"
#include <htslib/vcf.h>
#include <htslib/tbx.h>

}

