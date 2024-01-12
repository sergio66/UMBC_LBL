/*  
NAME

  hselect -- HITRAN record selection

SYNOPSIS

  1. Select a specified subset of a HITRAN database:

       hselect [-g id] [-r] [-i in] [-L v1] [-U v2] [-s n] <fin> [<fout>]

  2. Break a database into its componant gasses:

       hselect -b [-r] [-i in] [-L v1] [-U v2] [-s n] <fin> <dout>

OPTIONS

  -g id  select gas by HITRAN id (default is all gasses)
  -i in  select gas by isotope number (default is all isotopes)
  -L v1  lower bound for selected wavenumbers (default 0)
  -U v2  upper bound for selected wavenumbers (default 10000)
  -s n   select line strengths of at least n (default 0)
  -r     write output in a human-readable format (default is HITRAN)
  -b     break input database into separate gasses


DISCUSSION

  With no output spec, output is sent to stdout

  Multiple select options give the conjunction of the individual
  specifications.

  For the breakout option, -b, gasses in <dout> have filenames of
  the form g<n>.dat

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "hdefs.h"

#define MAXGAS 120

main (int argc, char *argv[]) {

  int i, j, k, c;
  struct HRAW hraw, hraw2;
  struct HSTR hstr, hstr2;
  struct HVAL hval, hval2;
  char *pname, *fin, *fout;
  char fgas[100];
  FILE *infid, *outfid;
  FILE *gasfid[MAXGAS];

  /* options defaults */
  int 	  gasid  = 0;
  int     isonum = 0;
  double  wL     = 0;
  double  wU     = 30000;
  double  smin   = 0;
  int     bflag  = 0;
  int     outfmt = 0;
  int     debug  = 1;

  extern int getopt();
  extern char *optarg;
  extern int   optind;
  extern char *sbase();

  pname = argv[0];

  while ((c=getopt(argc, argv, "g:i:L:U:s:brd")) != EOF)
    switch (c) {

    case 'g':    /* select a particular gas id */
      gasid = atoi(optarg);
      break;

    case 'i':    /* select a particular isotope */
      isonum = atoi(optarg);
      break;

    case 'L':    /* wavenumber lower bound */
      wL = atof(optarg);
      break;

    case 'U':    /* wavenumber upper bound */
      wU = atof(optarg);
      break;

    case 's':    /* specify a line strength threshold */
      smin = atof(optarg);
      break;

    case 'b':    /* break database into individual gasses */
      bflag = 1;
      break;

    case 'r':    /* readable output */
      outfmt = 1;
      break;

    case 'd':    /* debugging info */
      debug = 1;
      break;

    case '?':    /* bad command */
      exit(1);
      break;
    }

  if (argc - optind < 1) {
    fprintf(stderr, "%s: error -- no input file specified\n", pname);
    exit(1);
  }
  fin = argv[optind++];

  if (argc - optind == 1)
    fout = argv[optind++];
  else
    fout = "-";

  if (debug) {
    fprintf(stderr, "fin = %s, fout = %s\n", fin, fout);
    fprintf(stderr, "bflag = %d, gasid = %d, smin = %g, L = %g, U = %g\n", 
	    bflag, gasid, smin, wL, wU);
  }

  /* open input file
   */
  if (strcmp(fin, "-") == 0) 
    infid = stdin;
  else if ((infid = fopen(fin, "r")) == NULL) {
    fprintf(stderr, "%s: can't open input file %s", pname, fin);
    exit(0);
  }

  if (!bflag) {

    /*  SELECT MODE: write a selected subset of the input file 
     *  to the output file
     */
    if (strcmp(fout, "-") == 0) 
      outfid = stdout;
    else if ((outfid = fopen(fout, "w")) == NULL) {
      fprintf(stderr, "%s: can't open output file %s", pname, fout);
      exit(0);
    }

    while (fread(&hraw, sizeof(hraw), 1, infid) == 1) {

      raw2str(&hraw, &hstr);
      str2val(&hstr, &hval);

      if ((gasid == 0 || hval.IGAS == gasid) &&
	  (isonum == 0 || hval.ISO == isonum) &&
	  hval.WNUM >= wL &&
	  hval.WNUM <= wU &&
	  hval.STREN >= smin) {
	if (outfmt == 0) 
	  write_hraw(&hraw, outfid);
	else if (outfmt == 1)
	  write_hval(&hval, outfid);
      }
    }
    fclose(outfid);
  }
  else {

    /*  SPLIT MODE: separate input file into componant gasses, 
     *  one gas per file, in a common directory
     */
    if (strcmp(fout, "-") == 0) {
      fprintf(stderr, "%s: output directory can't be stdout!\n", pname);
      exit(127);
    }
    else if (mkdir(fout, 0755) == -1) {
      fprintf(stderr, 
	      "%s: can't create output directory %s\n", 
	      pname, fout);
      exit(127);
    }
    
     for (i=0; i<MAXGAS; i++) gasfid[i] = NULL;

     while (fread(&hraw, sizeof(hraw), 1, infid) == 1) {

       raw2str(&hraw, &hstr);
       str2val(&hstr, &hval);

       if (gasfid[hval.IGAS] == NULL) {

	 /* open output file for new gas  
	  */
	 sprintf(fgas, "%s/g%d.dat", fout, hval.IGAS);
          
	 if ((gasfid[hval.IGAS] = fopen(fgas, "w")) == NULL) {
	   fprintf(stderr, "%s: can't open output file %s\n", pname, fgas);
	   exit(0);
	 }
       }

       if ((isonum == 0 || hval.ISO == isonum) &&
	   hval.WNUM >= wL &&
	   hval.WNUM <= wU &&
	   hval.STREN >= smin) {
	 if (outfmt == 0) 
	   write_hraw(&hraw, gasfid[hval.IGAS]);
	 else if (outfmt == 1)
	   write_hval(&hval, gasfid[hval.IGAS]);
       }
     }
     for (i=0; i<MAXGAS; i++)
       if (gasfid[i] != NULL)
	 fclose(gasfid[i]);
  } 
}

