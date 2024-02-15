//  A. L. Delcher
//
//  File:  kmer-hash.hh
//
//  Last Modified:  5 Mar 2013
//
//  Declarations for a data structure to store kmers


#ifndef  __KMER_HASH_HH_INCLUDED
#define  __KMER_HASH_HH_INCLUDED

#include  "delcher.hh"


const int  DEFAULT_PREFIX_LEN = 10;
  // The number of prefix characters to use as a subscript to the list of
  // suffixes of kmers with that prefix
const int  MIN_KMER_LEN = 5;
  // The minimum kmer length allowed


struct Binary_Mer_t
{
  unsigned  prefix;
  __uint128_t  suffix;
};

template <class DT> struct Kmer_Hdr_t
{
  vector <__uint128_t>  suffix;
  vector <DT>  info;
};

template <class DT> class Kmer_Hash_t
{
private:
  DT * dummy;
  int  kmer_len;
  int  prefix_len, suffix_len;
  unsigned  prefix_first_mask, prefix_last_mask;
  __uint128_t  suffix_first_mask, suffix_last_mask;
  Kmer_Hdr_t <DT>  * hdr_list;

public:
  Kmer_Hash_t
    (int k, int p = DEFAULT_PREFIX_LEN);
  ~ Kmer_Hash_t ();
  void  Binary_To_Kmer
    (const Binary_Mer_t & bin, char kmer [], bool uppercase = false);
  void  Clear
    (void);
  int  Count_All
    (void);
  int  Count_All_Select
    (bool (* select) (const DT &));
  void  Dump_Kmers
    (FILE * fp);
  void  Dump_Kmers_Select
    (FILE * fp, void (* select) (FILE * p, const char * s, const DT &));
  DT *  Find
    (const char * kmer);
  DT *  Find
    (const Binary_Mer_t & bin);
  void  Find_All
    (const char * kmer, vector <DT> & list);
  void  Find_Or_Insert
    (const char * kmer, const DT & info);
  void  Find_Or_Insert
    (const Binary_Mer_t bin, const DT & info);
  void  Fwd_Shift_In
    (Binary_Mer_t & bin, char ch);
  void  Insert
    (const char * kmer, const DT & info);
  void  Insert
    (const Binary_Mer_t & bin, const DT & info);
  void  Kmer_To_Binary
    (const char * kmer, unsigned & prefix, __uint128_t & suffix);
  void  Kmer_To_Binary
    (const char * kmer, Binary_Mer_t & bin)
  {
    Kmer_To_Binary (kmer, bin . prefix, bin . suffix);
  }
  void  Rev_Shift_In
    (Binary_Mer_t & bin, char ch);
  void  Set_First_Ch
    (Binary_Mer_t & bin, char ch);
  void  Set_Last_Ch
    (Binary_Mer_t & bin, char ch);
};


unsigned  Char_To_Bits
  (char ch);


template <class DT> Kmer_Hash_t <DT> :: Kmer_Hash_t
  (int k, int p)

// Construct a kmer hash for kmers of length k, using an array to index
// to the position for entries with prefixes of length p.  Allocate
// space for s entries for each prefix

{
  long unsigned  n;
  int  max_suffix;

  kmer_len = k;
  if (k - 1 < p)
    p = k;
  prefix_len = p;
  suffix_len = kmer_len - prefix_len;

  // Using 2 bits per character check if suffixes will fit in
  // __uint128_t
  max_suffix = 4 * sizeof (__uint128_t);
  if (max_suffix < suffix_len)
    {
      sprintf (Clean_Exit_Msg_Line,
               "ERROR:  Kmer suffix length %d is too long; maximum is %d",
               suffix_len, max_suffix);
      Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
    }

  // Make masks to select bits for specific characters
  prefix_last_mask = 3;
    // last (tail) character of prefix
  prefix_first_mask = (prefix_last_mask << 2 * (prefix_len - 1));
    // first (front) character of prefix
  suffix_last_mask = 3;
    // last (tail) character of suffix
  suffix_first_mask = (suffix_last_mask << 2 * (suffix_len - 1));
    // first (front) character of suffix

  // Allocate an array of 2^(2p) entries, one for each possible prefix
  n = 1 << (2 * p);
  hdr_list = new Kmer_Hdr_t <DT> [n];
}


template <class DT> Kmer_Hash_t <DT> :: ~ Kmer_Hash_t ()

// Destroy this kmer hash

{
  delete [] hdr_list;
}


template <class DT> void  Kmer_Hash_t <DT> :: Binary_To_Kmer
  (const Binary_Mer_t & bin, char kmer [], bool uppercase)

// Convert binary representatin of kmer in bin to string format in kmer
// kmer must already have memory allocated to hold the result

{
  const char  ALPHA [] = "ACGT";
  const char  alpha [] = "acgt";
  unsigned  a, pref;
  __uint128_t  b, suff;
  int  i;

  pref = bin . prefix;
  for (i = 0; i < prefix_len; i ++)
    {
      a = pref & prefix_last_mask;
      kmer [prefix_len - i - 1] = (uppercase ? ALPHA [a] : alpha [a]);
      pref >>= 2;
    }

  suff = bin . suffix;
  for (i = 0; i < suffix_len; i ++)
    {
      b = suff & suffix_last_mask;
      kmer [kmer_len - i - 1] = (uppercase ? ALPHA [unsigned (b)] : alpha [unsigned (b)]);
      suff >>= 2;
    }

  kmer [kmer_len] = '\0';
}


template <class DT> void  Kmer_Hash_t <DT> :: Clear
  (void)

// Remove all entries from this hash

{
  int  i, n;

  n = 1 << (2 * prefix_len);
  for (i = 0; i < n; i ++)
    {
      hdr_list [i] . suffix . clear ();
      hdr_list [i] . info . clear ();
    }
}


template <class DT> int  Kmer_Hash_t <DT> :: Count_All
  (void)

// Return the number of kmers in this hash

{
  int  total = 0;
  int  i, n;

  n = 1 << (2 * prefix_len);
  for (i = 0; i < n; i ++)
    total += hdr_list [i] . suffix . size ();

  return total;
}


template <class DT> int  Kmer_Hash_t <DT> :: Count_All_Select
  (bool (* select) (const DT &))

// Return the number of kmers in this hash for which select is true

{
  int  total = 0;
  int  i, j, m, n;

  n = 1 << (2 * prefix_len);
  for (i = 0; i < n; i ++)
    {
      m = hdr_list [i] . info . size ();
      for (j = 0; j < m; j ++)
        if ((* select) (hdr_list [i] . info [j]))
          total ++;
    }

  return total;
}


template <class DT> void  Kmer_Hash_t <DT> :: Dump_Kmers
  (FILE * fp)

// Dump all kmers in this hash to fp

{
  Binary_Mer_t  bin;
  char  * s;
  unsigned  i, n;
  int  j, m;

  s = (char *) malloc (1 + kmer_len);

  n = 1 << (2 * prefix_len);
  for (i = 0; i < n; i ++)
    {
      bin . prefix = i;
      m = hdr_list [i] . suffix . size ();
      for (j = 0; j < m; j ++)
        {
          bin . suffix = hdr_list [i] . suffix [j];
          Binary_To_Kmer (bin, s);
          fprintf (fp, "%s\n", s);
        }
    }

  free (s);
}


template <class DT> void  Kmer_Hash_t <DT> :: Dump_Kmers_Select
  (FILE * fp, void (* select) (FILE * p, const char * s, const DT &))

// Run function select on all kmers in this hash, passing fp as the
// FILE * parameter of select.  Presumably, select chooses which kmers
// to print something about and what to print

{
  Binary_Mer_t  bin;
  char  * s;
  unsigned  i, n;
  int  j, m;

  s = (char *) malloc (1 + kmer_len);

  n = 1 << (2 * prefix_len);
  for (i = 0; i < n; i ++)
    {
      bin . prefix = i;
      m = hdr_list [i] . suffix . size ();
      for (j = 0; j < m; j ++)
        {
          bin . suffix = hdr_list [i] . suffix [j];
          Binary_To_Kmer (bin, s);
          (* select) (fp, s, hdr_list [i] . info [j]);
        }
    }

  free (s);
}


template <class DT> DT *  Kmer_Hash_t <DT> :: Find
  (const char * kmer)

// Find kmer in this hash and return a pointer to its
// corresponding info.  Return NULL if it's not found

{
  unsigned  prefix;
  __uint128_t  suffix;
  int  i, n;

  if (int (strlen (kmer)) != kmer_len)
    fprintf (stderr, "Oops:  strlen(kmer)= %d  kmer_len= %d\n",
             int (strlen (kmer)), kmer_len);
  assert (int (strlen (kmer)) == kmer_len);
  Kmer_To_Binary (kmer, prefix, suffix);

  vector <__uint128_t> &  sl = hdr_list [prefix] . suffix;

  n = sl . size ();
  for (i = 0; i < n; i ++)
    if (sl [i] == suffix)
      return & (hdr_list [prefix] . info [i]);

  return NULL;
}


template <class DT> DT *  Kmer_Hash_t <DT> :: Find
  (const Binary_Mer_t & bin)

// Find binary kmer bin in this hash and return a pointer to its
// corresponding info.  Return NULL if it's not found

{
  int  i, n;

  vector <__uint128_t> &  sl = hdr_list [bin . prefix] . suffix;

  n = sl . size ();
  for (i = 0; i < n; i ++)
    if (sl [i] == bin . suffix)
      return & (hdr_list [bin . prefix] . info [i]);

  return NULL;
}


template <class DT> void  Kmer_Hash_t <DT> :: Find_All
  (const char * kmer, vector <DT> & list)

// Find all matches of kmer in this hash and put all their
// corresponding info entries into list

{
  unsigned  prefix;
  __uint128_t  suffix;
  int  i, n;

  list . clear ();

  assert (int (strlen (kmer)) == kmer_len);
  Kmer_To_Binary (kmer, prefix, suffix);

  vector <__uint128_t> &  sl = hdr_list [prefix] . suffix;

  n = sl . size ();
  for (i = 0; i < n; i ++)
    if (sl [i] == suffix)
      list . push_back (hdr_list [prefix] . info [i]);

  return;
}


template <class DT> void  Kmer_Hash_t <DT> :: Find_Or_Insert
  (const char * kmer, const DT & info)

// Look for kmer in this hash.  If found, do nothing; otherwise,
// insert it with its associate info.

{
  Binary_Mer_t  bin;

  assert (int (strlen (kmer)) == kmer_len);
  Kmer_To_Binary (kmer, bin);

  Find_Or_Insert (bin, info);
}


template <class DT> void  Kmer_Hash_t <DT> :: Find_Or_Insert
  (const Binary_Mer_t bin, const DT & info)

// Search this hash for the kmer in bin.  If found, do nothing; otherwise,
// insert the kmer with its associate info.
  
{
  int  i, n;

  vector <__uint128_t> &  sl = hdr_list [bin . prefix] . suffix;

  n = sl . size ();
  for (i = 0; i < n; i ++)
    if (sl [i] == bin . suffix)
      return;

  sl . push_back (bin . suffix);
  hdr_list [bin . prefix] . info . push_back (info);
}


template <class DT> void  Kmer_Hash_t <DT> :: Fwd_Shift_In
  (Binary_Mer_t & bin, char ch)

// Remove the front (first) character from the kmer in bin;
// shift the other characters one position (2 bits) left;
// then add ch onto the end

{
  __uint128_t  a;

  bin . prefix &= ~ prefix_first_mask;
  bin . prefix <<= 2;
  a = bin . suffix & suffix_first_mask;
  a >>= (2 * (suffix_len - 1));
  bin . prefix |= unsigned (a);

  bin . suffix &= ~ suffix_first_mask;
  bin . suffix <<= 2;
  bin . suffix |= __uint128_t (Char_To_Bits (ch));
}


template <class DT> void  Kmer_Hash_t <DT> :: Insert
  (const char * kmer, const DT & info)

// Add kmer to this kmer hash with corresponding information in info
// Assume it's not already present

{
  __uint128_t  suff;
  unsigned  pref;
  printf("*****Printing strlen %d and kmer length %d",strlen(kmer), kmer_len);
  assert (int (strlen (kmer)) == kmer_len);
  Kmer_To_Binary (kmer, pref, suff);

  hdr_list [pref] . suffix . push_back (suff);
  hdr_list [pref] . info . push_back (info);
}


template <class DT> void  Kmer_Hash_t <DT> :: Insert
  (const Binary_Mer_t & bin, const DT & info)

// Add the kmer represented by bin to this kmer hash with
// corresponding information in info.  Do not check for duplicates

{
  hdr_list [bin . prefix] . suffix . push_back (bin . suffix);
  hdr_list [bin . prefix] . info . push_back (info);
}


template <class DT> void  Kmer_Hash_t <DT> :: Kmer_To_Binary
  (const char * kmer, unsigned & prefix, __uint128_t & suffix)

// Convert string kmer into a binary prefix (representing the first
// prefix_len characters) and a binary suffix (representing the rest).

{
  char  ch;
  int  i;

  prefix = suffix = 0;

  for (i = 0; i < prefix_len; i ++)
    {
      prefix <<= 2;
      ch = kmer [i];
      prefix |= Char_To_Bits (ch);
    }

  for ( ; i < kmer_len; i ++)
    {
      suffix <<= 2;
      ch = kmer [i];
      suffix |= (__uint128_t) (Char_To_Bits (ch));
    }
}


template <class DT> void  Kmer_Hash_t <DT> :: Rev_Shift_In
  (Binary_Mer_t & bin, char ch)

// Remove the tail (last) character from the kmer in bin;
// shift the other characters one position (2 bits) right;
// then add ch onto the front

{
  unsigned  a;
  __uint128_t  b;

  a = 3u & bin . prefix;
  b = __uint128_t (a);
  bin . suffix >>= 2;
  bin . suffix |= (b << 2 * (suffix_len - 1));

  bin . prefix >>= 2;
  a = Char_To_Bits (ch);
  bin . prefix |= (a << 2 * (prefix_len - 1));
}


template <class DT> void  Kmer_Hash_t <DT> :: Set_First_Ch
  (Binary_Mer_t & bin, char ch)

// Change the first (front) character of the kmer in bin to ch

{
  unsigned  a;

  a = Char_To_Bits (ch);

  bin . prefix &= (~ prefix_first_mask);
  bin . prefix |= (a << 2 * (prefix_len - 1));
}


template <class DT> void  Kmer_Hash_t <DT> :: Set_Last_Ch
  (Binary_Mer_t & bin, char ch)

// Change the last (tail) character of the kmer in bin to ch

{
  __uint128_t  a;

  a = Char_To_Bits (ch);

  bin . suffix &= (~ suffix_last_mask);
  bin . suffix |= a;
}



// ~~~ Below are from delcher.hh ~~~

#if 0
const int  MAX_ERROR_MSG_LEN = 1000;
  // Length of longest possible error message


extern char  Clean_Exit_Msg_Line [MAX_ERROR_MSG_LEN];
  // String to write error messages to before exiting
extern int  Verbose;
  // Flag to determine level of debugging output


const char *  Commatize
    (long int  n);
void  Clean_Exit
    (const char * msg, const char * src_fname = NULL, size_t line_num = 0);
FILE *  File_Open
    (const string & fname, const string & mode, const char * src_fname = NULL,
     size_t line_num = 0);
char  First_Non_Blank
    (const char * s);
int  Int_Power
    (int a, int b);
void  Make_Lower_Case
    (char * s);
void  Make_Upper_Case
    (char * s);
const char *  Num_Or_Max
    (int n, int mx = INT_MAX);
double  Percent
    (double a, double b);
const char *  Printable
    (bool b);
const char *  Printable
    (char * s);
double  Pseudo_Normal
    (void);
double  Ratio
    (double a, double b);
void  Reverse_String
    (char * s);
void  Reverse_String
    (string & s);
void *  Safe_calloc
    (size_t n, size_t len, const char * src_fname = NULL,
     size_t line_num = 0);
void *  Safe_malloc
    (size_t len, const char * src_fname = NULL, size_t line_num = 0);
void *  Safe_realloc
    (void * q, size_t len, const char * src_fname = NULL,
     size_t line_num = 0);
char *  Strip_Trailing
    (char * s, char ch);


template <class DT>  void  Incr_Limited
    (DT & A, DT limit);
template <class DT>  DT  Max
    (DT, DT);
template <class DT>  DT  Min
    (DT, DT);
template <class DT>  void  Swap
    (DT &, DT &);



template <class DT>  void  Incr_Limited
    (DT & A, DT limit)

// Increment  A  by 1, but only if it's less than  limit .

  {
   if  (A < limit)
       A ++;

   return;
  }



template <class DT>  DT  Max
    (DT A, DT B)

// Return the larger of  A  and  B .

  {
   if  (A > B)
       return  A;
     else
       return  B;
  }



template <class DT>  DT  Min
    (DT A, DT B)

// Return the smaller of  A  and  B .

  {
   if  (A < B)
       return  A;
     else
       return  B;
  }



template <class DT>  void  Reverse
    (vector <DT> & v)

// Reverse the order of entries in  v .

  {
   DT  s;
   int  i, j, n;

   n = v . size ();
   for  (i = 0, j = n - 1;  i < j;  i ++, j --)
     {
      s = v [i];
      v [i] = v [j];
      v [j] = s;
     }

   return;
  }



template <class DT>  void  Swap
    (DT & A, DT & B)

// Swap the values in  A  and  B .

  {
   DT  Save;

   Save = A;
   A = B;
   B = Save;

   return;
  }
#endif


#endif
