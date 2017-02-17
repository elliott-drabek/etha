//  A. L. Delcher
//
//  File:  kmer-repair.hh
//
//  Last Modified:  9 May 2013
//
//  Declarations for  kmer-repair.cc



#ifndef  __KMER_REPAIR_HH_INCLUDED
#define  __KMER_REPAIR_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "kmer-hash.hh"


const int  DEFAULT_HASH_PREFIX_CHARS = 10;
  // Default prefix length of kmers directly indexed in array
const int  DEFAULT_MAX_GAP = 500;
  // Default length of longest gap that will try to correct
const unsigned  DEFAULT_MAX_PATHS = 500;
  // Default maximum number of sequences that can be used for
  // extension or correction
const double  MAX_EXTENSION_ERATE = 0.12;
  // Maximum allowed alignment error rate in doing extensions
const int  MAX_EXTENSION_STEPS = 500;
  // Maximum number of steps to try in extending kmer matches
  // before the first match or after the last match
const int  MAX_LINE = 1000;
  // Max length of input lines/strings

static const char  COMPLEMENT_TABLE []
  = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
    " nnnnnnnnn*nn-.nnnnnnnnnnnnnnnnn"
    "nTVGHNNCDNNMNKNNNNYSANBWNRNnnnn_"
    "ntvghnncdnnmnknnnnysanbwnrnnnnnn";
static const char  FWD_ACGT [] = "ACGT";
static const char  REV_ACGT [] = "TGCA";


const int  KMER_INFO_FREQ_BITS = 14;

struct Correction_t
{
  int  lo, hi;
    // positions on ref sequence to be replaced, in 0-based gap coordinates
    // e.g., the first 3 bases are lo=0, hi=3; an insertion between characters
    // 3 and 4 (counting from 1) is lo=3, hi=3.
  string  replacement;
    // the string with which to replace that part of the ref sequence
};

struct Kmer_Info_t
{
  unsigned  freq : KMER_INFO_FREQ_BITS;
  unsigned  no_left_extension : 1;
  unsigned  no_right_extension : 1;

  Kmer_Info_t ()
  {
    freq = 0;
    no_left_extension = no_right_extension = 0;
  }
};

struct Kmer_Ref_t
{
  Binary_Mer_t  bin;
  int  len;
};


static int  Choose_Best_Left_Extension
  (const char * seq, int hi, vector <string> & path, int & dist,
   int & seq_lo, int & path_lo);
static int  Choose_Best_Path
  (const char * seq, int from, int to, vector <string> & path, int & dist);
static int  Choose_Best_Right_Extension
  (const char * seq, vector <string> & path, int & dist,
   int & seq_hi, int & path_hi);
char  Complement
  (char ch);
static int  Edit_Dist
  (const char * a, const char * b);
static bool  Find_Left_Extension
  (char * curr_path, int step, int max_steps, const char * fmer,
   const char * rmer, vector <string> & path, Kmer_Hash_t <Kmer_Info_t> * kmer_hash);
Kmer_Info_t  * Find_Mer
  (Kmer_Hash_t <Kmer_Info_t> * kmer_hash, const char * fmer, const char * rmer);
static bool  Find_Path
  (char * curr_path, int step, int max_steps, char * end_mer,
   vector <string> & path, Kmer_Hash_t <Kmer_Info_t> * kmer_hash);
static bool  Find_Right_Extension
  (char * curr_path, int step, int max_steps, const char * fmer,
   const char * rmer, vector <string> & path, Kmer_Hash_t <Kmer_Info_t> * kmer_hash);
static void  Get_Left_Extensions
  (const char * fmer, const char * rmer, vector <string> & path,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash, int max_extend);
static void  Get_Paths
  (const char * seq, int from, int to, vector <string> & path,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash);
static void  Get_Right_Extensions
  (const char * seq, int from, int seq_len, vector <string> & path,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash, int max_extend);
static void  Parse_Command_Line
  (int argc, char * argv []);
static bool  Path_Exists
  (const char * seq, int from, int to, int & len,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash);
static int  Prefix_Edit_Dist
  (const char * a, int m, const char * b, int n, int & a_hi, int & b_hi);
void  Reverse_Complement
  (char * s);
static int  Suffix_Edit_Dist
  (const char * a, int a_hi, const char * b, int & a_lo, int & b_lo);
static void  Usage
  (void);


#endif
