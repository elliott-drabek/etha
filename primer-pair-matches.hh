//  A. L. Delcher
//
//  File:  uncovered.hh
//
//  Last Modified:  10 Mar 2014
//
//  Declarations for  primer-pair-matches.cc



#ifndef  __PRIMER_PAIR_MATCHES_HH_INCLUDED
#define  __PRIMER_PAIR_MATCHES_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"


const int  MAX_LINE = 1000;
  // Length of longest input line

static const char  COMPLEMENT_TABLE []
  = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
    " nnnnnnnnn*nn-.nnnnnnnnnnnnnnnnn"
    "nTVGHNNCDNNMNKNNNNYSANBWNRNnnnn_"
    "ntvghnncdnnmnknnnnysanbwnrnnnnnn";


struct  Match_t
  {
    char  primer_dir;  // 'f' for fwd primer; 'r' for rev
    int  pos;  // position of match, counting from 1
    int  which;  // subscript of primer
  };


static char  Complement
  (char ch);
static void  Parse_Command_Line
  (int argc, char * argv []);
static void  Reverse_Complement
  (string & s);
static void  Usage
  (void);

#endif
