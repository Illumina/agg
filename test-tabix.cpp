#include "agg.h"
#include "utils.h"

extern "C" {
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/tbx.h"
}

// ichr,ifrom,ito are 0-based;
// returns -1 on error, 0 if the line is a comment line, 1 on success
static int _regions_parse_line(char *line, int ichr,int ifrom,int ito, char **chr,char **chr_end,int *from,int *to)
{
    *chr_end = NULL;

    if ( line[0]=='#' ) return 0;

    int k,l;    // index of the start and end column of the tab-delimited file
    if ( ifrom <= ito )
        k = ifrom, l = ito;
    else
        l = ifrom, k = ito;

    int i;
    char *se = line, *ss = NULL; // start and end
    char *tmp;
    for (i=0; i<=k && *se; i++)
    {
        ss = i==0 ? se++ : ++se;
        while (*se && *se!='\t') se++;
    }
    if ( i<=k ) return -1;
    if ( k==l )
    {
        *from = *to = strtol(ss, &tmp, 10);
        if ( tmp==ss ) return -1;
    }
    else
    {
        if ( k==ifrom )
            *from = strtol(ss, &tmp, 10);
        else
            *to = strtol(ss, &tmp, 10);
        if ( ss==tmp ) return -1;

        for (i=k; i<l && *se; i++)
        {
            ss = ++se;
            while (*se && *se!='\t') se++;
        }
        if ( i<l ) return -1;
        if ( k==ifrom )
            *to = strtol(ss, &tmp, 10);
        else
            *from = strtol(ss, &tmp, 10);
        if ( ss==tmp ) return -1;
    }

    ss = se = line;
    for (i=0; i<=ichr && *se; i++)
    {
        if ( i>0 ) ss = ++se;
        while (*se && *se!='\t') se++;
    }
    if ( i<=ichr ) return -1;
    *chr_end = se;
    *chr = ss;
    return 1;
}

int main(int argc,char **argv) {
    tbx_t *tbx;
    BGZF *fp;
    kstring_t s;
    int i;
    if ((tbx = tbx_index_load(argv[1])) == 0) die("problem loading index");
    if ((fp = bgzf_open(argv[1], "r")) == 0) die("problem opneing bed");
    int ichr=0, ifrom=1, ito=2,from,to,ret;
    char *chr, *chr_end;
    s.s = 0; s.l = s.m = 0;
    for (i =2; i < argc; ++i) {
        hts_itr_t *itr;
        if ((itr = tbx_itr_querys(tbx, argv[i])) == 0) continue;
        while (tbx_bgzf_itr_next(fp, tbx, itr, &s) >= 0) {
            puts(s.s);
            ret = _regions_parse_line(s.s, ichr,ifrom,abs(ito), &chr,&chr_end,&from,&to);
            cout << endl;
            cout<<" "<<from<<" "<<to<<endl;
        }

        tbx_itr_destroy(itr);
    }
    free(s.s);
    bgzf_close(fp);
    tbx_destroy(tbx);
}

