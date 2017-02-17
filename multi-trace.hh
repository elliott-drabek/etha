//  A. L. Delcher
//
//  File:  multi-trace.hh
//
//  Last Modified:  3 May 2013
//
//  Declarations for  multi-trace.cc



#ifndef  __MULTI_TRACE_HH_INCLUDED
#define  __MULTI_TRACE_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "kmer-hash.hh"


const int  DEFAULT_HASH_PREFIX_CHARS = 10;
  // Default prefix length of kmers directly indexed in array
const int  MAX_LINE = 1000;
  // Max length of input lines/strings

static const char  COMPLEMENT_TABLE []
  = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
    " nnnnnnnnn*nn-.nnnnnnnnnnnnnnnnn"
    "nTVGHNNCDNNMNKNNNNYSANBWNRNnnnn_"
    "ntvghnncdnnmnknnnnysanbwnrnnnnnn";
static const char  FWD_ACGT [] = "ACGT";
static const char  REV_ACGT [] = "TGCA";


const int  KMER_INFO_FREQ_BITS = 15;

struct Kmer_Info_t
{
  unsigned  freq : KMER_INFO_FREQ_BITS;
  unsigned  hit : 1;

  Kmer_Info_t ()
  {
    freq = hit = 0;
  }
};


char  Complement
  (char ch);
Kmer_Info_t  * Find_Mer
  (Kmer_Hash_t <Kmer_Info_t> * kmer_hash, const char * fmer, const char * rmer);
bool  Has_Positive_Hit
  (const Kmer_Info_t & info);
static void  Parse_Command_Line
  (int argc, char * argv []);
void  Reverse_Complement
  (char * s);
static void  Usage
  (void);


#endif
