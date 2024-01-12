/* 

% NAME
%   read_hitran - read HITRAN data, mex version
%
% SYNOPSIS
%   s = read_hitran(v1, v2, str, gid, hsrc);
%
% INPUTS
%   v1   - wavenumber lower bound
%   v2   - wavenumber upper bound
%   str  - line strength lower bound
%   gid  - HITRAN gas id 
%   hsrc - HITRAN data file or directory
%
% OUTPUTS
%   s.igas    - molecule ID number
%   s.iso     - isotope number
%   s.wnum    - line wavenumber
%   s.stren   - line strength
% 
%   s.tprob   - transition probability
%   s.abroad  - air-broad half width
%   s.sbroad  - self-broad half width
%   s.els     - lower-state energy
%   s.abcoef  - coeff of temp dep of ABROAD
% 
%   s.tsp     - transition shift due to pressure 
%   s.iusgq   - upper state global quanta index
%   s.ilsgq   - lower state global quanta index
%   s.uslq    - upper state local quanta
%   s.bslq    - lower state local quanta
%   s.ai      - accuracy indices
%   s.ref     - indices for lookup of references
%   s.flag    - flag
%   s.swus    - statistical weight, upper state
%   s.swls    - statistical weight, lower state
% 
% NOTES
%   The optional input parameter hsrc can be either a HITRAN data file,
%   or a directory of HITRAN files organized gas types.  For the latter
%   case, the individual files are assumed to have have names g<n>.dat,
%   for gas <n>.  It is possible (but very slow) to specify the whole,
%   unsplit HITRAN database file as hsrc, and select gasses by gid; the
%   usual way to use the procedure is to specify a directory for hsrc,
%   and select gasses by gasid.  
% 
%   The default name for hsrc is "hitran.dat"; this is convenient for
%   setting the gas directory with a symlink.
%
% AUTHOR
%   H. Motteler, 15 Dec 06

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>

#include "mex.h"
#include "matrix.h"

#include "hdefs.h"


/* read_hitran(wL, wU, smin, gasid, fin)
 * 
 * reads a file of HITRAN format (160 char) records, selects
 * records, and return selected data as a linked list of HVAL
 * structures
 *  
 */

int nrec;

struct HVAL *read_hitran (double wL, double wU, 
			   double smin, int gasid, char* fin) 
{
  struct HRAW hraw;
  struct HSTR hstr;
  struct HVAL hval; 
  struct HVAL *hvalP0 = NULL;
  struct HVAL *hvalP1 = NULL;
  struct HVAL *hvalP2 = NULL;
  char cbuf[120];
  FILE *infid;
  DIR *dtst;

  /* option defaults */
  int isonum = 0;

  /* build input filename if we were passed a directory
   */ 
  if ((dtst = opendir(fin)) != NULL) {
    closedir(dtst);
    sprintf(cbuf, "%s/g%d.dat", fin, gasid);
    fin = cbuf;
  }

  /* open input file
   */
  if ((infid = fopen(fin, "r")) == NULL)
    mexErrMsgTxt("can't open input file");

  /* loop on hitran records
   */
  nrec = 0;
  while (fread(&hraw, sizeof(hraw), 1, infid) == 1) {

    raw2str(&hraw, &hstr);
    str2val(&hstr, &hval);

    /* assume data sorted in wavenumber order */
    if (hval.WNUM > wU) break; 

    if ((gasid == 0 || hval.IGAS == gasid) &&
	(isonum == 0 || hval.ISO == isonum) &&
	hval.WNUM >= wL &&
	hval.WNUM <= wU &&
	hval.STREN >= smin) {

      hvalP2 = hvalP1;
      hvalP1 = (struct HVAL *) mxCalloc(1, sizeof(struct HVAL));
      if (hvalP1 == NULL)
	mexErrMsgTxt("mxCalloc failed");

      if (hvalP2)
	hvalP2->next = hvalP1;

      *hvalP1 = hval;
      hvalP1->next = NULL;

      if (hvalP0 == NULL)
	hvalP0 = hvalP1;

      nrec++;
    }
  }
  fclose(infid);
  return hvalP0;
}


/* source is char string, dest is matlab char string
 */
void mxncpy(mxChar *a, char *b, int n) {
  int i=0;
  mxChar *p;
  char *q;
  for (p=a, q=b; i < n; i++) *p++ = *q++;
}


/* mex gateway procedure
 *
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  /* input args */
  double wL, wU;
  double smin, gasr;
  char* fbuf;

  /* output args */
  mxArray *outstr;
  mxArray *outmxa[19];
  double  *outptr[19];
  const char *fnames[] = {
    "igas",   /*  0 */
    "iso",    /*  1 */
    "wnum",   /*  2 */
    "stren",  /*  3 */
    "tprob",  /*  4 */
    "abroad", /*  5 */
    "sbroad", /*  6 */
    "els",    /*  7 */
    "abcoef", /*  8 */
    "tsp",    /*  9 */
    "iusgq",  /* 10 */
    "ilsgq",  /* 11 */
    "uslq",   /* 12 */
    "bslq",   /* 13 */
    "ai",     /* 14 */
    "ref",    /* 15 */
    "flag",   /* 16 */
    "swus",   /* 17 */
    "swls"};  /* 18 */

  /* character field lengths */
  int cflen[] = { 15, 15, 15, 15, 6, 12, 1 };

  /* msc variables */
  int ndims[2];
  int i, j, m, n, flen, gasi;
  struct HVAL * hval;   
  mxArray *mxp1[1], *mxp2[1];

  /* Check for proper number of arguments 
   */
  if (nrhs != 4 && nrhs != 5)
    mexErrMsgTxt("read_hitran() requires 4 or 5 input arguments");
  else if (nlhs > 1)
    mexErrMsgTxt("read_hitran(): too many output arguments");
  
  wL = mxGetScalar(prhs[0]);
  wU = mxGetScalar(prhs[1]);
  smin = mxGetScalar(prhs[2]);
  gasr = mxGetScalar(prhs[3]);
  gasi = gasr;

  /* filename string */
  if (nrhs == 5) {
    m = mxGetM(prhs[4]);
    n = mxGetN(prhs[4]);
    flen = m * n + 1;
    fbuf = (char *) mxCalloc(flen, sizeof(char));
    mxGetString(prhs[4], fbuf, flen);
  }
  else
    fbuf = "hitran.dat";


  /* call read_hitran() to do the real work */
  hval = read_hitran(wL, wU, smin, gasi, fbuf);

  /* create real arrays for output */
  for (i=0; i<10; i++) {
    outmxa[i] = mxCreateDoubleMatrix(nrec, 1, mxREAL);
    if (outmxa[i] == NULL)
      mexErrMsgTxt("can't create mx output array");
  }

  /* creat char arrays for output */
  for (i=10; i<17; i++) {
    ndims[0] = cflen[i-10];
    ndims[1] = nrec;
    outmxa[i] = mxCreateCharArray(2, ndims);    
    if (outmxa[i] == NULL)
      mexErrMsgTxt("can't create mx output array");
  }

  /* last two fields are real */
  for (i=17; i<19; i++) {
    outmxa[i] = mxCreateDoubleMatrix(nrec, 1, mxREAL);
    if (outmxa[i] == NULL)
      mexErrMsgTxt("can't create mx output array");
  }

  /* get pointers to actual data-space in the mxArrays */
  for (i=0; i<19; i++) 
    outptr[i] = mxGetPr(outmxa[i]);

  /* copy HVAL records to matlab output arrays */
  for (i=0; i<nrec; i++) {

    *(outptr[0] + i)  =  hval->IGAS;
    *(outptr[1] + i)  =  hval->ISO;
    *(outptr[2] + i)  =  hval->WNUM;
    *(outptr[3] + i)  =  hval->STREN;
    *(outptr[4] + i)  =  hval->TPROB;
    *(outptr[5] + i)  =  hval->ABROAD;
    *(outptr[6] + i)  =  hval->SBROAD;
    *(outptr[7] + i)  =  hval->ELS;
    *(outptr[8] + i)  =  hval->ABCOEF;
    *(outptr[9] + i)  =  hval->TSP;

    mxncpy((mxChar *) outptr[10] + i*15, hval->IUSGQ, 15);
    mxncpy((mxChar *) outptr[11] + i*15, hval->ILSGQ, 15);
    mxncpy((mxChar *) outptr[12] + i*15, hval->USLQ, 15);
    mxncpy((mxChar *) outptr[13] + i*15, hval->BSLQ, 15);
    mxncpy((mxChar *) outptr[14] + i*6,  hval->AI,  6);
    mxncpy((mxChar *) outptr[15] + i*12, hval->REF, 12);
    mxncpy((mxChar *) outptr[16] + i*1, hval->FLAG, 1);

    *(outptr[17] + i) = hval->SWUS;
    *(outptr[18] + i) = hval->SWLS;

    hval = hval->next;
  }    

  /* transpose character arrays */
  for (i=10; i<17; i++) {
    mxp1[0] = outmxa[i];
    mexCallMATLAB(1, mxp2, 1, mxp1, "transpose");
    outmxa[i] = mxp2[0];
  }

  /* get structure to bundle results */
  outstr = mxCreateStructMatrix(1, 1, 19, fnames);

  for (i=0; i<19; i++) 
    mxSetFieldByNumber(outstr, 0, i, outmxa[i]);

   plhs[0] = outstr;
}

