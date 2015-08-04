#include <iostream>
#include <stdlib.h>     /* atoi */
using namespace std;

int main ()
{
    string chrom;
    int mingq = 20;
    int a=-1,b,dp,gq;
    string tmp_gq,tmp_b,tmp_dp;
    int block_start=-1,block_end,block_gq,block_dp;
    string prev_chrom;
    while(cin>>chrom && cin>>a && cin>>tmp_b && cin>>tmp_gq && cin>>tmp_dp) {
        cout <<"\t"<<chrom <<"\t"<<a <<"\t"<<b <<"\t"<<tmp_gq <<"\t"<<dp<<endl;//debugging
        if(tmp_gq==".") gq=0;
        else gq = atoi(tmp_gq.c_str());
        if(tmp_b==".") b=a;
        else b = atoi(tmp_b.c_str());

        if(block_start>-1) {
            if( (a-1)==block_end && gq>=mingq ) {//extend current block.
                block_end=b;
                if(block_gq>gq) block_gq=gq;
                if(block_dp>dp) block_dp=dp;
            }
            else {//print block and reset
                cout << prev_chrom <<"\t"<< block_start <<"\t"<< block_end <<"\t"<< block_gq <<"\t"<< block_dp<<endl;
                block_start=-1;	
            }
        }      
        if(block_start==-1 && gq>=mingq) {//intialise
            block_start=a;
            block_end=b;
            block_gq=gq;
            block_dp=dp;
        }
        prev_chrom=chrom;
    }

    return 0;
}
