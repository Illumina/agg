#include "agg.h"
#include "kseq.h"
#include "kstring.h"
#include <deque>
#include <set>
#include <string>         // std::string
#include <locale>         // std::locale, std::toupper
#include <algorithm>
#include <iostream>
#include "needle.h"
#include "hts_utils.h"

KSTREAM_INIT(gzFile, gzread, 16384)


extern "C" {
#include "htslib/hts.h"
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "vcfnorm.h"
}


#define ERR_REF_MISMATCH    -1
#define CHECK_REF_WARN 1
