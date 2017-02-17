//  A. L. Delcher
//
//  File:  unitig.hh
//
//  Last Modified:  14 Oct 2013
//
//  Declarations for  unitig.cc



#ifndef  __UNITIG_HH_INCLUDED
#define  __UNITIG_HH_INCLUDED


#include  "delcher.hh"


const int  MAX_LINE = 1000;
  // Max length of input lines/strings

struct Olap_t
{
  int  a_lo, a_hi, b_lo, b_hi;
  int  a_sub, b_sub;
};

struct Seq_t
{
  char * id;
  int  len : 30;
  unsigned  contained : 1;
  unsigned  used : 1;
  vector <Olap_t>  left_olap, right_olap;

  Seq_t ()  // constructor
  {
    contained = used = 0;
  }

  bool Has_Right_Olap_To (int p)
  {
    int  i, n;

    n = right_olap . size ();
    for (i = 0; i < n; i ++)
      if (right_olap [i] . b_sub == p)
        return true;
    return false;
  }

  bool Has_Left_Olap_To (int p)
  {
    int  i, n;

    n = left_olap . size ();
    for (i = 0; i < n; i ++)
      if (left_olap [i] . b_sub == p)
        return true;
    return false;
  }
};


static int  Find_Or_Insert
  (const char * b_id, int b_len, vector <Seq_t> & seq);
static void  Find_Unitigs
  (vector <Seq_t> & seq);
static void  Parse_Command_Line
  (int argc, char * argv []);
static void  Remove_Contained_Olaps
  (vector <Seq_t> & seq);
static void  Remove_Transitive_Olaps
  (vector <Seq_t> & seq);
static void  Usage
  (void);


#endif
