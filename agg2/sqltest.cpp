#include <stdio.h>
#include <sqlite3.h>
#include "utils.h"

extern "C" {
#include "htslib/synced_bcf_reader.h"
}


static int callback(void *NotUsed, int argc, char **argv, char **azColName){
    int i;
    for(i=0; i<argc; i++){
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return 0;
}

int createDB(sqlite3 *db) {
    char *zErrMsg = 0;
    if(SQLITE_OK!=sqlite3_exec(db, "CREATE TABLE variants (chromosome ENUM, position UNSIGNED INT, ref TEXT, alt TEXT, ac UNSIGNED INT, npass SINGLE, qual SINGLE );", callback, 0, &zErrMsg)) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }

    if(SQLITE_OK!=sqlite3_exec(db,"CREATE INDEX IF NOT EXISTS pos_idx ON variants(chromosome, position);" ,callback, NULL, &zErrMsg)) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }

    return(0);
}


int updateDB(sqlite3 *db,char *input_vcf) {
    bcf_srs_t *reader= bcf_sr_init();
    if(!bcf_sr_add_reader (reader, input_vcf)) die("problem opening "+((string)input_vcf));
    bcf_hdr_t *header = reader->readers[0].header;  
    bcf1_t *line;
    char *zErrMsg = 0;
    sqlite3_stmt *stmt;  
    int rc,ac;
    int ngt,mgt_arr = 0, *gt_arr = NULL;

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
    sqlite3_prepare_v2(db, "INSERT INTO variants VALUES (?1,?2,?3,?4,?5,?6,?7);", -1, &stmt, NULL); 
    int n=0;
    while(bcf_sr_next_line (reader)) {  
        line = bcf_sr_get_line(reader,0);    
        ngt = bcf_get_genotypes(header, line, &gt_arr, &mgt_arr);
        int a = bcf_gt_allele(gt_arr[0]);
        int b = bcf_gt_allele(gt_arr[1]);  
        sqlite3_bind_text(stmt, 1, (char *)bcf_seqname(header,line), -1, SQLITE_STATIC);//chrom
        sqlite3_bind_int(stmt, 2, line->pos);//position
        sqlite3_bind_text(stmt, 3,line->d.allele[0], -1, SQLITE_STATIC);//REF		    
        sqlite3_bind_int(stmt, 6, bcf_has_filter(header, line, ".") );//pass?
        sqlite3_bind_double(stmt, 7, line->qual);

        for(unsigned int i=1;i<line->n_allele;i++) {       
            int ac = (a==i) + (b==i);
            sqlite3_bind_text(stmt, 4,line->d.allele[i], -1, SQLITE_STATIC);//ALT
            sqlite3_bind_int(stmt, 5, ac);//AC
            rc = sqlite3_step(stmt); 
            if (rc != SQLITE_DONE) {
                die("ERROR inserting data: %s\n"+(string)sqlite3_errmsg(db));
            }
        }

        sqlite3_clear_bindings(stmt);
        sqlite3_reset(stmt);
        n++;
    }

    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);
    bcf_sr_destroy(reader);
    return(n);
}


int main(int argc, char **argv){
    sqlite3 *db;
    int rc;

    if( argc!=3 ){
        fprintf(stderr, "Usage: %s DATABASE SQL-STATEMENT\n", argv[0]);
        return(1);
    }
    rc = sqlite3_open(argv[1], &db);
    if( rc ){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return(1);
    }
    createDB(db);
    updateDB(db,argv[2]);

    sqlite3_close(db);
    return 0;
}
