//  A. L. Delcher
//
//  File:  multi-walk.hh
//
//  Last Modified:  5 Apr 2013
//
//  Declarations for  multi-walk.cc



#ifndef  __MULTI_WALK_HH_INCLUDED
#define  __MULTI_WALK_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "kmer-hash.hh"


const int  DEFAULT_NUM_STEPS = 200;
  // Default number of steps to take in walk
const int  MAX_LINE = 1000;
  // Max length of input lines/strings
const int  MIN_OLAP = 5;
  // A match overlapping another by this much or more on both strings
  // is discarded

static const char  COMPLEMENT_TABLE []
  = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
    " nnnnnnnnn*nn-.nnnnnnnnnnnnnnnnn"
    "nTVGHNNCDNNMNKNNNNYSANBWNRNnnnn_"
    "ntvghnncdnnmnknnnnysanbwnrnnnnnn";
static const char  FWD_ACGT [] = "ACGT";
static const char  REV_ACGT [] = "TGCA";


const int  KMER_INFO_FREQ_BITS = 14;

struct Kmer_Info_t
{
  unsigned  freq : KMER_INFO_FREQ_BITS;
  unsigned  used_fwd : 1;    // true if this version of kmer is used
  unsigned  used_rc : 1;     // true if reverse complement of kmer is used

  Kmer_Info_t ()
  {
    freq = 0;
    used_fwd = used_rc = 0;
  }
};

struct Result_Hash_Entry_t
{
  unsigned  empty : 1;
  int  pos : 31;  // start position of kmer in string
};

struct Stack_Entry_t
{
  int  pos, freq;
  int  r_lo, r_hi;  // give active portion of result_string
  int  frame;
  char  dir;
  string  fmer, rmer;
  unsigned char  open_frame;
    // Use 1 bit for each frame: 0x1 for frame 0; 0x2 for frame 1; 0x4 for frame 2
    // so, in general, (1 << f) for frame f
};

struct Adj_t
{
  int  id;
  int  mult;  // multiplicity of the edge to this node
  Adj_t  * next;
};

struct Node_t
{
  int  id;
  int  in_degr, out_degr;
    // total of multiplicities of edges into and out of this node
  Adj_t  * to, * from;
    // first in list of nodes this node has connections to and from, respectively
  bool  active;
};


static bool  Check_Stops
  (vector <char *> & stop_pattern, char * fmer, char dir, int & pattern_sub);
char  Complement
  (char ch);
static bool  Is_Stop
  (const char * codon);
static int  Num_Open
  (unsigned char open_frame);
static void  Parse_Command_Line
  (int argc, char * argv []);
void  Print_If_Used
  (FILE * fp, const char * s, const Kmer_Info_t & info);
static void  Print_Result
  (FILE * fp, char * s, int lo, int hi, char * id, unsigned ver, const char * expl);
void  Reverse_Complement
  (char * s);
static int  Set_Open_Frame
  (unsigned char open_frame, const char * s, int offset);
static void  Usage
  (void);


#endif
