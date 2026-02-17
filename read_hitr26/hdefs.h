

/* "raw" concatenated input format
 */
struct HRAW {
  char IGAS[2];
  char ISO[1];
  char WNUM[12];
  char STREN[10];
  char TPROB[10];
  char ABROAD[5];
  char SBROAD[5];
  char ELS[10];
  char ABCOEF[4];
  char TSP[8];
  char IUSGQ[15];
  char ILSGQ[15];
  char USLQ[15];
  char BSLQ[15];
  char AI[6];
  char REF[12];
  char FLAG[1];
  char SWUS[7];
  char SWLS[7];
  char NL[1];
};


/* string format with null separators
 * (intermediate step for conversions)
 */
struct HSTR {
  char IGAS[3];
  char ISO[2];
  char WNUM[13];
  char STREN[11];
  char TPROB[11];
  char ABROAD[6];
  char SBROAD[6];
  char ELS[11];
  char ABCOEF[5];
  char TSP[9];
  char IUSGQ[16];
  char ILSGQ[16];
  char USLQ[16];
  char BSLQ[16];
  char AI[7];
  char REF[13];
  char FLAG[2];
  char SWUS[8];
  char SWLS[8];
};


/* hitran values in native data types
 * (the last four fields are left as strings)
 */
struct HVAL{
  int    IGAS;
  int    ISO;
  double WNUM;
  double STREN;
  double TPROB;
  double ABROAD;
  double SBROAD;
  double ELS;
  double ABCOEF;
  double TSP;
  char   IUSGQ[16];
  char   ILSGQ[16];
  char   USLQ[16];
  char   BSLQ[16];
  char   AI[7];
  char   REF[13];
  char   FLAG[2];
  double SWUS;
  double SWLS;
  struct HVAL *next;
};


/* prototypes
 */
void raw2str (struct HRAW *hraw, struct HSTR *hstr);
void str2raw (struct HSTR *hstr, struct HRAW *hraw);
void str2val (struct HSTR *hstr, struct HVAL *hval);
void val2str (struct HVAL *hval, struct HSTR *hstr);
void write_hval (struct HVAL *hval, FILE *fid);
void write_hraw (struct HRAW *hraw, FILE *fid);
void write_hstr (struct HSTR *hstr, FILE *fid);

