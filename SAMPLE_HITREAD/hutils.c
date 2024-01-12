
/*  string utilities for HITRAN records
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "hdefs.h"


/*  copy exactly n characters from s2 to s1 and
 *  add a null at the end of s1, at position n+1.
 */
void scp(char *a, char *b, int n) {
  int i=0;
  char *p, *q;
  for (p=a, q=b; i < n; i++) *p++ = *q++;
  *p = 0;
}


/* copy a raw record to a null-terminated string record
 */
void raw2str (struct HRAW *hraw, struct HSTR *hstr) {

  scp(hstr->IGAS,   hraw->IGAS,    2);
  scp(hstr->ISO,    hraw->ISO,     1);
  scp(hstr->WNUM,   hraw->WNUM,    12);
  scp(hstr->STREN,  hraw->STREN,   10);
  scp(hstr->TPROB,  hraw->TPROB,   10);
  scp(hstr->ABROAD, hraw->ABROAD,  5);
  scp(hstr->SBROAD, hraw->SBROAD,  5);
  scp(hstr->ELS,    hraw->ELS,     10);
  scp(hstr->ABCOEF, hraw->ABCOEF,  4);
  scp(hstr->TSP,    hraw->TSP,     8);

  scp(hstr->IUSGQ,  hraw->IUSGQ,   15);
  scp(hstr->ILSGQ,  hraw->ILSGQ,   15);
  scp(hstr->USLQ,   hraw->USLQ,    15);
  scp(hstr->BSLQ,   hraw->BSLQ,    15);
  scp(hstr->AI,     hraw->AI,      6);
  scp(hstr->REF,    hraw->REF,     12);

  scp(hstr->FLAG,   hraw->FLAG,    1);
  scp(hstr->SWUS,   hraw->SWUS,    7);
  scp(hstr->SWLS,   hraw->SWLS,    7);
} 


/* copy a null-terminated string record back to a raw record
 */
void str2raw (struct HSTR *hstr, struct HRAW *hraw) {

  scp(hraw->IGAS,   hstr->IGAS,    2);
  scp(hraw->ISO,    hstr->ISO,     1);
  scp(hraw->WNUM,   hstr->WNUM,    12);
  scp(hraw->STREN,  hstr->STREN,   10);
  scp(hraw->TPROB,  hstr->TPROB,   10);
  scp(hraw->ABROAD, hstr->ABROAD,  5);
  scp(hraw->SBROAD, hstr->SBROAD,  5);
  scp(hraw->ELS,    hstr->ELS,     10);
  scp(hraw->ABCOEF, hstr->ABCOEF,  4);
  scp(hraw->TSP,    hstr->TSP,     8);

  scp(hraw->IUSGQ,  hstr->IUSGQ,   15);
  scp(hraw->ILSGQ,  hstr->ILSGQ,   15);
  scp(hraw->USLQ,   hstr->USLQ,    15);
  scp(hraw->BSLQ,   hstr->BSLQ,    15);
  scp(hraw->AI,     hstr->AI,      6);
  scp(hraw->REF,    hstr->REF,     12);

  scp(hraw->FLAG,   hstr->FLAG,    1);
  scp(hraw->SWUS,   hstr->SWUS,    7);
  scp(hraw->SWLS,   hstr->SWLS,    7);

  hraw->NL[0] = '\n';
} 


/* copy a null-terminated string record to a value record
 */
void str2val (struct HSTR *hstr, struct HVAL *hval) {

  hval->IGAS    =  atoi(hstr->IGAS);
  hval->ISO     =  atoi(hstr->ISO);
  hval->WNUM    =  atof(hstr->WNUM);
  hval->STREN   =  atof(hstr->STREN);
  hval->TPROB   =  atof(hstr->TPROB);
  hval->ABROAD  =  atof(hstr->ABROAD);
  hval->SBROAD  =  atof(hstr->SBROAD);
  hval->ELS     =  atof(hstr->ELS);
  hval->ABCOEF  =  atof(hstr->ABCOEF);
  hval->TSP     =  atof(hstr->TSP);

  scp(hval->IUSGQ, hstr->IUSGQ, 15);
  scp(hval->ILSGQ, hstr->ILSGQ, 15);
  scp(hval->USLQ,  hstr->USLQ, 15);
  scp(hval->BSLQ,  hstr->BSLQ, 15);
  scp(hval->AI,    hstr->AI,  6);
  scp(hval->REF,   hstr->REF, 12);

  scp(hval->FLAG,  hstr->FLAG, 1);
  hval->SWUS = atof(hstr->SWUS);
  hval->SWLS = atof(hstr->SWLS);
}


/* copy a value record to a null-terminated string record
 */
static void verr(char *s) {
/*     fprintf(stderr, "val2str(): bad %s length\n", s); */
/*     exit(127); */
}

void val2str (struct HVAL *hval, struct HSTR *hstr) {

char b[8];

  if (sprintf(hstr->IGAS,   "%2d",    hval->IGAS)   != 2)  verr("IGAS");
  if (sprintf(hstr->ISO,    "%1d",    hval->ISO)    != 1)  verr("ISO");
  if (sprintf(hstr->WNUM,   "%12.6f", hval->WNUM)   != 12) verr("WNUM");
  if (sprintf(hstr->STREN,  "%10.3e", hval->STREN)  != 10) verr("STREN");
  if (sprintf(hstr->TPROB,  "%10.3e", hval->TPROB)  != 10) verr("TPROB");

  /* a hack to strip off the leading zero for a fortran-style F5.4 */

  if (sprintf(b, "%6.4f",  hval->ABROAD) != 6)  verr("ABROAD");
  scp(hstr->ABROAD, b+1, 5);

  if (sprintf(b, "%6.4f",  hval->SBROAD) != 6)  verr("SBROAD");
  scp(hstr->SBROAD, b+1, 5);

  if (sprintf(hstr->ELS,    "%10.4f", hval->ELS)    != 10) verr("ELS");
  if (sprintf(hstr->ABCOEF, "%4.2f",  hval->ABCOEF) != 4)  verr("ABCOEF");
  if (sprintf(hstr->TSP,    "%8.6f",  hval->TSP)    != 8)  verr("TSP");

  if (sprintf(hstr->IUSGQ,  "%15s",   hval->IUSGQ)  != 15)  verr("IUSGQ");
  if (sprintf(hstr->ILSGQ,  "%15s",   hval->ILSGQ)  != 15)  verr("ILSGQ");
  if (sprintf(hstr->USLQ,   "%15s",   hval->USLQ)   != 15)  verr("USLQ");
  if (sprintf(hstr->BSLQ,   "%15s",   hval->BSLQ)   != 15)  verr("BSLQ");
  if (sprintf(hstr->AI,     "%6s",    hval->AI)     !=  6)  verr("AI");
  if (sprintf(hstr->REF,    "%12s",   hval->REF)    != 12)  verr("REF");

  if (sprintf(hstr->FLAG,   "%1s",    hval->FLAG)   != 1)  verr("FLAG");
  if (sprintf(hstr->SWUS,   "%7.1f",  hval->SWUS)   != 7)  verr("SWUS");
  if (sprintf(hstr->SWLS,   "%7.1f",  hval->SWLS)   != 7)  verr("SWLS");
}


/*  write a "value record" in readable format
 */
void write_hval (struct HVAL *hval, FILE *fid) {

  fprintf(fid, 
	  "%2d %1d %12.6f %10.3e %10.3e %6.4f %6.4f %10.4f %4.2f %8.6f \
	   %15s %15s %15s %15s %6s %12s %1s %f7.1 %f7.1\n",
	 hval->IGAS,
	 hval->ISO,
	 hval->WNUM,
	 hval->STREN,
	 hval->TPROB,
	 hval->ABROAD,
	 hval->SBROAD,
	 hval->ELS,
	 hval->ABCOEF,
	 hval->TSP,
	 hval->IUSGQ,
	 hval->ILSGQ,
	 hval->USLQ,
	 hval->BSLQ,
	 hval->AI,
  	 hval->REF,
	 hval->FLAG,
	 hval->SWUS,
	 hval->SWLS
	 );
}


/* write a raw record
 */
void write_hraw (struct HRAW *hraw, FILE *fid) {

  fwrite(hraw, sizeof(struct HRAW), 1, fid);
}


/* write a string-format record 
 * (mainly for testing purposes...)
 */
void write_hstr (struct HSTR *hstr, FILE *fid) {
  int i;
  char c;
  for (i=0; i<sizeof(struct HSTR); i++) {
    c = ((char *) hstr)[i];
    if (c == 0) 
      putc(' ', fid);
    else
      putc(c, fid); 
  }
  putc('\n', fid);
}

