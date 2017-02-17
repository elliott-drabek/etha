//  A. L. Delcher
//
//  File:  kmer-repair.cc
//
//  Last Modified:  9 May 2013
//
//  This program reads a list of kmers with counts, in the format produced
//  by 'jellyfish dump -c'.  It then reads a DNA multi-fasta sequence file
//  and looks for matches of its kmers.  For regions between kmer matches,
//  the paths of kmers connecting the matches are compared to the sequence
//  and the closest one to the sequence is stitched into the sequence to
//  correct/repair it.  Might also repair off the ends of the extreme kmer
//  matches where there is a unique kmer path.


#include  "kmer-repair.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

int  Hash_Prefix_Chars = DEFAULT_HASH_PREFIX_CHARS;
  // Prefix length of kmers directly indexed in array
const char  * Kmer_Count_Filename;
  // Name of file from which to read kmer counts
int  Kmer_Len = -1;
  // Length of kmers
int  Max_Gap = DEFAULT_MAX_GAP;
  // Longest gap that will attempt to correct
unsigned  Max_Paths = DEFAULT_MAX_PATHS;
  // Maximum number of sequences that can be used for extension or
  // correction.  If more than this many sequences are found, no correction
  // is attempted.
const char  * Sequence_Filename;
  // Name of file with DNA fasta sequences


int  main
  (int argc, char * argv [])

{
  FILE  * sequence_fp, * kmer_count_fp;
  Kmer_Hash_t <Kmer_Info_t>  * kmer_hash;
  string  seq_string, seq_hdr;
  Kmer_Info_t  info;
  char  fmer [MAX_LINE], rmer [MAX_LINE], last_fmer [MAX_LINE];
  char  s [MAX_LINE];
  const char  * seq;
  unsigned  max_count;
  int  ct, seq_len, max_extend;
  int  i, j, n;

  Verbose = 0;

  Parse_Command_Line (argc, argv);

  if (strcmp (Sequence_Filename, "-") == 0)
    sequence_fp = stdin;
  else
    sequence_fp = File_Open (Sequence_Filename, "r");

  if (strcmp (Kmer_Count_Filename, "-") == 0)
    {
      kmer_count_fp = stdin;
      if (sequence_fp == stdin)
        {
          sprintf (Clean_Exit_Msg_Line,
                   "ERROR:  Sequence file and kmer count file cannot both be '-'");
          Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
        }
    }
  else
    kmer_count_fp = File_Open (Kmer_Count_Filename, "r");

  max_count = (1 << KMER_INFO_FREQ_BITS) - 1;
  printf ("max_count = %d\n", max_count);

  fprintf (stderr, "Reading kmer counts ...\n");
  i = 0;
  while (fscanf (kmer_count_fp, "%s %d", s, & ct) == 2)
    {
      if (i == 0)
        {
          Kmer_Len = strlen (s);
          // Create the kmer hash
          kmer_hash = (Kmer_Hash_t <Kmer_Info_t> *) new Kmer_Hash_t <Kmer_Info_t> (Kmer_Len, Hash_Prefix_Chars);
        }
      
      if (max_count < unsigned (ct))
        info . freq = max_count;
      else
        info . freq = ct;

      for (j = 0; j < Kmer_Len; j ++)
        s [j] = toupper (s [j]);
      kmer_hash -> Insert (s, info);

      i ++;
      if (i % 10000000 == 0)
        fprintf (stderr, " %4d million counts read\n", i / 1000000);

#if 0
      if (72000000 < i)
        {
          printf ("Truncated reading kmer counts after %d lines\n", i);
          break;
        }
#endif
    }
  fprintf (stderr, "Total counts read = %d\n", i);

  if (kmer_count_fp != stdin)
    fclose (kmer_count_fp);

  // Read and process the fasta sequences
  while (Fasta_Read (sequence_fp, seq_string, seq_hdr))
    {
      vector <string>  path;
      vector <Correction_t>  correct;
      Correction_t  corr;
      Kmer_Info_t  * p;
      double  erate;
      int  dist, last_match;

      seq_len = seq_string . length ();
      if (seq_len < Kmer_Len)
        {
          fprintf (stderr,
                   "Skipping sequence %s with length %d shorter than kmer length %d\n",
                   seq_hdr . c_str (), seq_len, Kmer_Len);
          continue;
        }

      for (i = 0; i < seq_len; i ++)
        seq_string [i] = toupper (seq_string [i]);
      seq = seq_string . c_str ();

      printf ("\n#%s  len=%d\n", seq_hdr . c_str (), seq_len);

      strncpy (fmer, seq, Kmer_Len);
      fmer [Kmer_Len] = '\0';
      strcpy (rmer, fmer);
      Reverse_Complement (rmer);
      
      last_match = -1;
      for (i = 1; i <= seq_len - Kmer_Len + 1; i ++)
        {
          char  tag [MAX_LINE];

          // Don't allow matches to N's.  They can happen because N's
          // are converted to A's in the kmer hash.
          if (strchr (fmer, 'N') == NULL)
            p = Find_Mer (kmer_hash, fmer, rmer);
          else
            p = NULL;

          if (p == NULL)
            ct = 0;
          else
            ct = p -> freq;

          if (p == NULL)
            sprintf (tag, "unmatched");
          else
            sprintf (tag, "ct= %d", ct);
          printf ("%5d  %s  %s\n", i, fmer, tag);

          if (p != NULL)
            {
              if (last_match == -1)
                { // this is the first match and there are prior unmatched kmers
                  // see if there are extensions left from it that can be used
                  // to repair the sequence
                  if (1 < i)
                    {
                      max_extend = int (1.1 * (i - 1) + 10);
                        // have i-1 unmatched kmers; allow extra characters
                        // in extension in case of indels
                      
                      Get_Left_Extensions (fmer, rmer, path, kmer_hash, max_extend);
                      n = path . size ();
                    }
                  else
                    n = 0;

                  if (Max_Paths < unsigned (n))
                    {
                      printf ("  Found more than %d extensions--ignoring\n", Max_Paths);
                      path . clear ();
                      n = 0;
                    }
                  else if (1 < i)
                    printf ("  Found %d left extensions\n", n);

                  if (0 < n)
                    {
                      bool  made_correction = false;
                      int  seq_lo, path_lo;

                      j = Choose_Best_Left_Extension (seq, i - 2, path,
                                                      dist, seq_lo, path_lo);
                      if (j < 0)
                        erate = DBL_MAX;
                      else
                        erate = (1.0 * dist) / (path [j] . length () - path_lo);
                      if (path_lo == 0 && erate <= MAX_EXTENSION_ERATE)
                        {
                          corr . lo = seq_lo;
                          corr . hi = i - 1;
                          corr . replacement = path [j] . substr (path_lo,
                                                 path [j] . length () - path_lo);
                          correct . push_back (corr);
                          made_correction = true;
                        }
                      if (0 < Verbose)
                        printf ("Best left extension is path j= %d  dist= %d"
                                "  erate= %.2f%%\n",
                                j, dist, 100.0 * erate);
                      if (0 <= j && 0 < Verbose)
                        {
                          string  s, t;

                          t = seq;
                          s = t . substr (seq_lo, i - 1 - seq_lo);
                          printf ("Replace:  %s\n", s . c_str ());
                          s = path [j] . substr (path_lo,
                                                 path [j] . length () - path_lo);
                          printf ("   with:  %s\n", s . c_str ());
                          if (made_correction)
                            printf ("  Correction made\n");
                          else
                            printf ("  Correction NOT made\n");
                        }
                    }
                }
              else if (i != last_match + 1)
                {
                  int  len;
                  int  gap = i - last_match - 1;

                  printf ("  Unmatched gap of length %d\n", gap);
                  printf ("  Search paths from %s\n", last_fmer);
                  printf ("                 to %s\n", fmer);
                  
                  if (gap <= Max_Gap)
                    {
                      // First check if a path exists so we don't waste
                      // too much time enumerating hopeless paths
                      // i values start at one; string positions start at zero
//**ALD  Changed here
                      if (Path_Exists (seq, last_match - 1, i - 1, len, kmer_hash))
                        {
                          Get_Paths (seq, last_match - 1, i - 1, path, kmer_hash);
                          n = path . size ();
                          if (Max_Paths < unsigned (n))
                            {
                              printf ("  Found more than %d paths--ignoring\n", Max_Paths);
                              path . clear ();
                              n = 0;
                            }
                          else
                            printf ("  Found %d paths\n", n);
                        }
                      else
                        {
                          path . clear ();
                          n = 0;
                        }

                      if (0 < n)
                        {
                          int  path_sep, read_sep;

                          j = Choose_Best_Path (seq, last_match - 1, i - 1, path, dist);
                          printf ("best path= %d  score= %d\n", j, dist);
                          corr . lo = last_match + Kmer_Len - 1;
                          corr . hi = i - 1;
                          read_sep = corr . hi - corr . lo;
                            // offset between start & end kmers on read to be corrected
                          path_sep = path [j] . length () - 2 * Kmer_Len;
                            // offset between those kmers on correct path
                          if (0 <= read_sep && 0 <= path_sep)
                            { // simple replacement
                              corr . replacement
                                = path [j] . substr (Kmer_Len, path_sep);
                            }
                          else if (read_sep < path_sep)
                            { // read_sep is negative and need insertion
                              corr . hi = corr . lo;
                              corr . replacement
                                = path [j] . substr (Kmer_Len, path_sep - read_sep);
                            }
                          else if (path_sep < read_sep)
                            { // path_sep is negative and need deletion
                              corr . hi = corr . lo + read_sep - path_sep;
                              corr . replacement = "";
                            }
                          else
                            { // this can't happen--the read & path sequences would be
                              // the same and there wouldn't be a gap to fix
                              fprintf (stderr, "read_sep= %d  path_sep= %d\n",
                                       read_sep, path_sep);
                              assert (false);
                            }
                          
                          if (0 < Verbose)
                            {
                              char  * old;


                              old = (char *) malloc (1 + corr . hi - corr . lo);
                              strncpy (old, seq + corr . lo, corr . hi - corr . lo);
                              old [corr . hi - corr . lo] = '\0';
                              printf ("Replace pos %d..%d\n"
                                      "     >%s<\n"
                                      "with >%s<\n",
                                      corr . lo, corr . hi, old,
                                      corr . replacement . c_str ());
                              free (old);
                            }
                          correct . push_back (corr);
                        }
                    }
                  else
                    printf ("  Gap too long--skipping\n");
                }
              last_match = i;
              strcpy (last_fmer, fmer);
            }

          // Advance to next character
          memmove (fmer, fmer + 1, Kmer_Len - 1);
          memmove (rmer + 1, rmer, Kmer_Len - 1);
          fmer [Kmer_Len - 1] = seq [i + Kmer_Len - 1];
          rmer [0] = Complement (fmer [Kmer_Len - 1]);
        }

      if (0 < last_match && last_match < seq_len - Kmer_Len + 1)
        {
          j = seq_len - Kmer_Len + 1 - last_match;
          max_extend = int (1.1 * j + 10);
            // have j unmatched kmers; allow extra characters
            // in extension in case of indels

          printf ("Try right extension  last_match= %d  seq_len= %d  max_extend= %d\n",
                  last_match, seq_len, max_extend);

          Get_Right_Extensions (seq, last_match - 1, seq_len, path, kmer_hash,
                                max_extend);

          n = path . size ();
          if (Max_Paths < unsigned (n))
            {
              printf ("  Found more than %d extensions--ignoring\n", Max_Paths);
              path . clear ();
              n = 0;
            }
          else
            printf ("  Found %d right extensions\n", n);

          if (0 < n)
            {
              bool  made_correction = false;
              int  seq_hi, path_hi;

              j = Choose_Best_Right_Extension (seq + last_match + Kmer_Len - 1,
                                               path, dist, seq_hi, path_hi);

              printf ("j= %d\n", j);
              if (j < 0)
                erate = DBL_MAX;
              else
                erate = (2.0 * dist) / (seq_hi + path_hi);
              if (erate <= MAX_EXTENSION_ERATE)
                {
                  corr . lo = last_match + Kmer_Len - 1;
                  corr . hi = corr . lo + seq_hi;
                  corr . replacement = path [j] . substr (0, path_hi);
                  correct . push_back (corr);
                  made_correction = true;
                }
              if (0 < j && 0 < Verbose)
                {
                  string  s, t;

                  printf ("Best right extension is path j= %d  dist= %d"
                          "  erate= %.2f%%\n",
                          j, dist, 100.0 * erate);
                  t = seq;
                  s = t . substr (last_match + Kmer_Len - 1, seq_hi);
                  printf ("Replace:  %s\n", s . c_str ());
                  s = path [j] . substr (0, path_hi);
                  printf ("   with:  %s\n", corr . replacement . c_str ());
                  if (made_correction)
                    printf ("  Correction made\n");
                  else
                    printf ("  Correction NOT made\n");
                }
            }
        }

      if (0 < Verbose)
        printf ("Number of correction segments is %lu\n", correct . size ());

      // Output corrected sequence
      string  corr_seq;
      int  prev;

      n = correct . size ();
      prev = 0;
      for (i = 0; i < n; i ++)
        {
          corr_seq . append (seq_string . substr (prev, correct [i] . lo - prev));
          corr_seq . append (correct [i] . replacement);
          prev = correct [i] . hi;
        }
      if (prev < seq_len)
        corr_seq . append (seq_string . substr (prev, seq_len - prev));
      if (0 < n)
        seq_hdr . append ("  corrected");
      Fasta_Print (stdout, corr_seq . c_str (), seq_hdr . c_str ());
    }

  if (sequence_fp != stdin)
    fclose (sequence_fp);

  return  0;
}


static int  Choose_Best_Left_Extension
  (const char * seq, int hi, vector <string> & path, int & dist,
   int & seq_lo, int & path_lo)

// Find suffix edit distances from strings in path to string seq [0 .. hi].
// Return the subscript of the one with the highest identity match.
// Set dist to the suffix edit distance of that string and set seq_lo and
// path_lo to the start positions of that match in seq and in the path
// string, respectively.

{
  double  best_erate, erate;
  int  best = -1;
  int  s_lo, p_lo;
  int  d, i, n;

  if (0 < Verbose)
    printf ("Choose_Best_Left_Extension:\n");
  dist = INT_MAX;  // impossibly high value
  best_erate = 2.0;  // impossibly high

  n = path . size ();
  for (i = 0; i < n; i ++)
    {
      d = Suffix_Edit_Dist (seq, hi, path [i] . c_str (), s_lo, p_lo);
      if (0 < 1 + hi - s_lo + path [i] . length () - p_lo)
        erate = (2.0 * d) / (1 + hi - s_lo + path [i] . length () - p_lo);
      else
        erate = 2.0;
      if (0 < Verbose)
        {
          printf ("path= %d  d= %d  erate= %6.4f\n", i, d, erate);
          printf ("  hi= %d  s_lo= %d  p_lo= %d\n", hi, s_lo, p_lo);
        }
      if (erate < best_erate)
        {
          best_erate = erate;
          dist = d;
          best = i;
          seq_lo = s_lo;
          path_lo = p_lo;
        }
    }

  return best;
}


static int  Choose_Best_Path
  (const char * seq, int from, int to, vector <string> & path, int & dist)

// Find the edit distance from each string in path to the substring of
// seq from position from to position to + Kmer_Len.  Set dist to the edit
// distance of that string and return its subscript.  In case of ties,
// return the first string found.

{
  char  * ref;
  int  best = -1, len;
  int  d, i, n;

  if (0 < Verbose)
    printf ("Choose_Best_Path:\n");
  dist = len = to + Kmer_Len - from;
  ref = (char *) malloc (len + 1);

  if (0 < Verbose)
    printf ("from= %d  to= %d  len= %d  seq_len= %lu\n", from, to, len, strlen (seq));
  strncpy (ref, seq + from, len);
  ref [len] = '\0';

  n = path . size ();
  for (i = 0; i < n; i ++)
    {
      d = Edit_Dist (ref, path [i] . c_str ());
      if (0 < Verbose)
        printf ("path= %d  d= %d\n", i, d);
      if (d < dist)
        {
          dist = d;
          best = i;
        }
    }

  free (ref);

  return best;
}


char  Complement
  (char ch)

// Returns the DNA complement of  ch

{
  return  COMPLEMENT_TABLE [unsigned (ch)];
}


static int  Choose_Best_Right_Extension
  (const char * seq, vector <string> & path, int & dist,
   int & seq_hi, int & path_hi)

// Find prefix edit distances from strings in path to string seq.
// Return the subscript of the one with the highest identity match.
// Set dist to the prefix edit distance of that string and set seq_hi
// and path_hi to the end positions of that match in seq and in the
// path string, respectively.

{
  double  best_erate, erate;
  int  best = -1;
  int  s_hi, p_hi, path_len, seq_len;
  int  d, i, n;

  if (0 < Verbose)
    printf ("Choose_Best_Right_Extension:\n");
  dist = INT_MAX;  // impossibly high value
  best_erate = 2.0;  // impossibly high

  seq_len = strlen (seq);
  n = path . size ();
  for (i = 0; i < n; i ++)
    {
      path_len = path [i] . length ();
      d = Prefix_Edit_Dist (seq, seq_len, path [i] . c_str (), path_len, s_hi, p_hi);
      if (0 < s_hi + p_hi)
        erate = (2.0 * d) / (s_hi + p_hi);
      else
        erate = 2.0;
      if (0 < Verbose)
        {
          printf ("path= %d  d= %d  erate= %6.4f\n", i, d, erate);
          printf ("  seq_len= %d  path_len= %d  s_hi= %d  p_hi= %d\n",
                  seq_len, path_len, s_hi, p_hi);
        }
      if ((s_hi == seq_len || p_hi == path_len) && erate < best_erate)
        {
          best_erate = erate;
          dist = d;
          best = i;
          seq_hi = s_hi;
          path_hi = p_hi;
        }
    }

  return best;
}


static int  Edit_Dist
  (const char * a, const char * b)

// Return the minimum number of single-character insertions, deletions and/or
// substitutions needed to convert string a into string b.  Uses the
// Vishkin-Schieber algorithm

{
  static short int  * space = NULL;
  static int  save_sz = -1;
  int  from, max_errs, sz;
  int  d, e, i, j, m, n;

  if (1 < Verbose)
    printf ("Edit_Dist:\n");
  m = strlen (a);
  n = strlen (b);

  max_errs = int (0.5 + 0.10 * Max (m, n));  // Allow up to 10% error
  sz = (1 + max_errs) * (1 + max_errs);

  if (save_sz < sz)
    {
      space = (short int *) SAFE_REALLOC (space, sz * sizeof (short int));
      save_sz = sz;
    }

  if (1 < Verbose)
    {
      printf ("m= %d  n= %d  max_errs= %d  sz= %d\n", m, n, max_errs, sz);
      Fasta_Print (stdout, a, "a-string");
      Fasta_Print (stdout, b, "b-string");
    }

#define  S(e,d)  space [(e) * ((e) + 1) + (d)]
  // Treat space as a pyramidal array with one cell on row 0, 3 on row 1,
  // 5 on row 2, ....  Cells on row e are indexed by d = -e ... +e

  for (i = 0; i < m && i < n && a [i] == b [i]; i ++)
    ;
  S (0, 0) = i;
  if (i == m && i == n)
    return  0;  // strings are identical

  if (1 < Verbose)
    printf ("e= %d  d= %d  S= %d\n", 0, 0, S (0, 0));
  for (e = 1; e <= max_errs; e ++)
    for (d = -e; d <= e; d ++)
      {
        i = -1;
        if (-e < d && d < e)
          {
            i = 1 + S (e - 1, d);
            from = d;
          }
        if (- e < d - 1 && i < (j = S (e - 1, d - 1)))
          {
            i = j;
            from = d - 1;
          }
        if (d + 1 < e && i < (j = 1 + S (e - 1, d + 1)))
          {
            i = j;
            from = d + 1;
          }
        while (i < m && i + d < n && a [i] == b [i + d])
          i ++;
        S (e, d) = i;
        if (i == m && i + d == n)
          return e;
        if (1 < Verbose)
          printf ("e= %d  d= %d  from= %d  S= %d  sub= %d\n",
                  e, d, from, S (e, d), e * (e + 1) + d);
      }

  return max_errs + 1;
}


static bool  Find_Left_Extension
  (char * curr_path, int step, int max_steps, const char * fmer,
   const char * rmer, vector <string> & path, Kmer_Hash_t <Kmer_Info_t> * kmer_hash)

// Recursively search for strings of overlapping kmers from kmer_hash
// that extend the string in (curr_path + step), going left from the
// the kmer pair fmer, rmer (forward and reverse strands).  Stop when
// step = 0.  If a non-empty path is found, it is appended to the
// vector of strings in path.  If the number of paths exceeds global
// Max_Paths, stop and return false; otherwise, return true.

{
  char  a [Kmer_Len + 1], b [Kmer_Len + 1];
  Kmer_Info_t  * p;
  bool  extended, ok;
  int  i;

  if (1 < Verbose)
    printf ("Find_Left_Extension:  step= %d  max_steps= %d\n  fmer= %s\n",
            step, max_steps, fmer);

  ok = true;
  extended = false;
  if (0 < step)
    {
      strncpy (a + 1, fmer, Kmer_Len - 1);
      a [Kmer_Len] = '\0';
      strncpy (b, rmer + 1, Kmer_Len - 1);
      b [Kmer_Len] = '\0';
  
      for (i = 0; i < 4 && ok; i ++)
        {
          a [0] = FWD_ACGT [i];
          b [Kmer_Len - 1] = REV_ACGT [i];
          p = Find_Mer (kmer_hash, a, b);
          if (p == NULL)
            continue;
          else
            {
              if (1 < Verbose)
                printf ("  step= %d  added_ch= %c  freq= %d\n", step, FWD_ACGT [i],
                        p -> freq);
              curr_path [step - 1] = FWD_ACGT [i];
              ok = Find_Left_Extension (curr_path, step - 1, max_steps, a, b, path,
                                        kmer_hash);
              extended = true;
            }
        }
    }

  if (! ok)
    return false;

  if (! extended && step != max_steps)
    { // Add curr_path to path
      path . push_back (curr_path + step);
      if (0 < Verbose)
        {
          printf ("Extension ended at step %d\n", step);
          printf ("Path extension is:\n  %s\n", curr_path + step);
        }
    }

  return (path . size () <= Max_Paths);
}

  
Kmer_Info_t  * Find_Mer
  (Kmer_Hash_t <Kmer_Info_t> * kmer_hash, const char * fmer, const char * rmer)

// Search for the kmer with forward-strand sequence fmer and reverse-strand sequence
// rmer in kmer_hash.  If found return a pointer to its corresponding information
// entry; otherwise, return NULL

{
  const char  * s;

  if (strcmp (fmer, rmer) < 0)
    s = fmer;
  else
    s = rmer;

  return kmer_hash -> Find (s);
}


static bool  Find_Path
  (char * curr_path, int step, int max_steps, char * end_mer,
   vector <string> & path, Kmer_Hash_t <Kmer_Info_t> * kmer_hash)

// Recursively search for a string of overlapping kmers from kmer_hash
// that extends the string curr_path until it reaches the kmer in
// end_mer.  step is the start position of the last kmer in curr_path.
// max_steps is the maximum number of steps to try before giving up.
// If a path is found, it is appended to the vector of strings in path.
// If the number of paths exceeds global Max_Paths, stop and return false;
// otherwise, return true.

{
  char  fmer [MAX_LINE], rmer [MAX_LINE];
  Kmer_Info_t  * p;
  bool  ok;
  int  i;

  if (1 < Verbose)
    printf ("Start step= %d  max_steps= %d\n", step, max_steps);
  if (step == max_steps)
    return true;

  curr_path [step + Kmer_Len] = '\0';
  strcpy (fmer, curr_path + step);

  if (1 < Verbose)
    printf ("fmer= %s\n", fmer);

  if (strcmp (fmer, end_mer) == 0)
    { // found the end_mer
      path . push_back (curr_path);
      if (0 < Verbose)
        {
          printf ("Found end_mer at step %d\n", step);
          printf ("Path extension is:\n  %s\n", curr_path + Kmer_Len);
        }
      return (path . size () <= Max_Paths);
    }

  // Look for extensions
  memmove (fmer, fmer + 1, Kmer_Len - 1);
  fmer [Kmer_Len - 1] = 'A';
  strcpy (rmer, fmer);
  Reverse_Complement (rmer);

  ok = true;
  for (i = 0; i < 4 && ok; i ++)
    {
      fmer [Kmer_Len - 1] = FWD_ACGT [i];
      rmer [0] = REV_ACGT [i];
      p = Find_Mer (kmer_hash, fmer, rmer);
      if (p == NULL)
        continue;
      if (1 < Verbose)
        printf ("  step= %d  added_ch= %c  freq= %d\n", step, FWD_ACGT [i],
                p -> freq);
      curr_path [step + Kmer_Len] = FWD_ACGT [i];
      ok = Find_Path (curr_path, step + 1, max_steps, end_mer, path, kmer_hash);
    }

  if (1 < Verbose)
    printf ("Back from step= %d\n", step);

  return ok;
}


static bool  Find_Right_Extension
  (char * curr_path, int step, int max_steps, const char * fmer,
   const char * rmer, vector <string> & path, Kmer_Hash_t <Kmer_Info_t> * kmer_hash)

// Recursively search for strings of overlapping kmers from kmer_hash
// that extend the string in curr_path (which has length step), going
// right from the the kmer pair fmer, rmer (forward and reverse
// strands) until no further extension is possible.  If a non-empty
// path is found, it is appended to the vector of strings in path.
// Stop and append the path if step = max_steps is reached.
// If the number of paths exceeds global Max_Paths, stop and return false;
// otherwise, return true.

{
  char  a [Kmer_Len + 1], b [Kmer_Len + 1];
  Kmer_Info_t  * p;
  bool  extended, ok;
  int  i;

  if (1 < Verbose)
    printf ("Find_Right_Extension:  step= %d  max_steps= %d\n  fmer= %s\n",
            step, max_steps, fmer);

  ok = true;
  extended = false;
  if (step < max_steps)
    {
      strncpy (a, fmer + 1, Kmer_Len - 1);
      a [Kmer_Len] = '\0';
      strncpy (b + 1, rmer, Kmer_Len - 1);
      b [Kmer_Len] = '\0';
  
      for (i = 0; i < 4 && ok; i ++)
        {
          a [Kmer_Len - 1] = FWD_ACGT [i];
          b [0] = REV_ACGT [i];
          p = Find_Mer (kmer_hash, a, b);
          if (p == NULL)
            continue;

          if (1 < Verbose)
            printf ("  step= %d  added_ch= %c  freq= %d\n", step, FWD_ACGT [i],
                    p -> freq);
          curr_path [step] = FWD_ACGT [i];
          ok = Find_Right_Extension (curr_path, step + 1, max_steps, a, b, path,
                                kmer_hash);
          extended = true;
        }
    }

  if (! ok)
    return false;

  if (! extended && 0 < step)
    { // Add curr_path to path
      curr_path [step] = '\0';
      path . push_back (curr_path);
      if (0 < Verbose)
        {
          printf ("Extension ended at step %d\n", step);
          printf ("Path extension is:\n  %s\n", curr_path);
        }
    }

  return (path . size () <= Max_Paths);
}


static void  Get_Left_Extensions
  (const char * fmer, const char * rmer, vector <string> & path,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash, int max_extend)

// Find strings made of consecutive kmers from kmer_hash going left
// from the kmer represented in fmer & rmer (forward and reverse
// strands of the kmer) until there is no next kmer.  Put the
// resulting strings in path.  max_extend is the maximum length string
// to produce.

{
  char  curr_path [MAX_EXTENSION_STEPS + 1];
  bool  ok;
  int  max_steps;

  max_steps = Min (MAX_EXTENSION_STEPS, max_extend);

  if (0 < Verbose)
    printf ("Get_Left_Extensions:  fmer= %s  max_steps= %d\n", fmer, max_steps);

  path . clear ();

  curr_path [max_steps] = '\0';

  ok = Find_Left_Extension (curr_path, max_steps, max_steps, fmer, rmer,
                       path, kmer_hash);

  if (0 < Verbose)
    {
      int  i, n = path . size ();

      printf ("Found %d left extensions  ok= %c\n", n, Bool_To_Char (ok));
      for (i = 0; i < n; i ++)
        printf ("  %s\n", path [i] . c_str ());
    }

  return;
}


static void  Get_Paths
  (const char * seq, int from, int to, vector <string> & path,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash)

// Find strings of consecutive kmers from kmer_hash going from the
// kmer starting at position from in seq to the kmer starting at
// position to.  Put the paths in path.

{
  char  * curr_path, * end_mer;
  bool  ok;
  int  max_steps;

  path . clear ();

  max_steps = int (10 + 1.1 * (to - from));
  curr_path = (char *) malloc (max_steps + Kmer_Len + 1);
  end_mer = (char *) malloc (Kmer_Len + 1);

  strncpy (curr_path, seq + from, Kmer_Len);
  curr_path [Kmer_Len] = '\0';
  strncpy (end_mer, seq + to, Kmer_Len);
  end_mer [Kmer_Len] = '\0';

  ok = Find_Path (curr_path, 0, max_steps, end_mer, path, kmer_hash);

  free (curr_path);
  free (end_mer);
  
  return;
}


static void  Get_Right_Extensions
  (const char * seq, int from, int seq_len, vector <string> & path,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash, int max_extend)

// Find strings of consecutive kmers from kmer_hash starting at
// the kmer beginning at seq [from] and going right until there
// is no next kmer.  Put the resulting strings in path.  max_extend is
// the maximum length string to produce.
  
{
  char  fmer [Kmer_Len + 1], rmer [Kmer_Len + 1];
  char  curr_path [MAX_EXTENSION_STEPS + 1];
  bool  ok;
  int  max_steps;

  max_steps = Min (MAX_EXTENSION_STEPS, max_extend);

  path . clear ();

  strncpy (fmer, seq + from, Kmer_Len);
  fmer [Kmer_Len] = '\0';
  strcpy (rmer, fmer);
  Reverse_Complement (rmer);

  if (0 < Verbose)
    printf ("Get_Right_Extensions:  from= %d\n  fmer= %s  max_steps= %d\n",
            from, fmer, max_steps);

  curr_path [0] = '\0';

  ok = Find_Right_Extension (curr_path, 0, max_steps, fmer, rmer, path, kmer_hash);
  
  if (0 < Verbose)
    {
      int  i, n = path . size ();

      printf ("Found %d right extensions\n", n);
      for (i = 0; i < n; i ++)
        printf ("  %s\n", path [i] . c_str ());
    }

  return;
}


static void  Parse_Command_Line
  (int argc, char * argv [])

// Get options and parameters from command line with  argc
// arguments in  argv [0 .. (argc - 1)] .

{
  bool  errflg = false;
  int  ch;

  optarg = NULL;

  while  (! errflg
          && ((ch = getopt (argc, argv, "g:hV:x:")) != EOF))
    switch  (ch)
      {
      case  'g' :
        Max_Gap = strtol (optarg, NULL, 10);
        break;

      case  'h' :
        Usage ();
        exit (EXIT_SUCCESS);

      case  'V' :
        Verbose = strtol (optarg, NULL, 10);
        break;

      case  'x' :
        Max_Paths = strtol (optarg, NULL, 10);
        if (Max_Paths < 1)
          {
            fprintf (stderr, "ERROR:  Bad -x option (max paths) value %d\n",
                     Max_Paths);
            errflg = true;
          }
        break;

      default :
        errflg = true;
      }

  if  (errflg || optind > argc - 2)
    {
      Usage ();
      exit (EXIT_FAILURE);
    }

  Sequence_Filename = argv [optind ++];
  Kmer_Count_Filename = argv [optind ++];

  return;
}


static bool  Path_Exists
  (const char * seq, int from, int to, int & len,
   Kmer_Hash_t <Kmer_Info_t> * kmer_hash)

// Find the shortest path of kmers in kmer_hash going from the
// kmer starting at position from in seq to the kmer starting at
// position to.  If there is a path, return true and set len to
// the length of the shortest path; otherwise, return false
// and set len to -1.  Note that we don't check if the start kmer
// (at position from) or end kmer (at position to) are in kmer_hash.

{
  Kmer_Hash_t <Kmer_Info_t>  used (Kmer_Len, Hash_Prefix_Chars);
  std :: queue <Kmer_Ref_t>  kmer_queue;
  Kmer_Info_t  dummy, * p;
  Kmer_Ref_t  kref;
  char  * start_mer, * end_mer, * fmer, * rmer;
  bool  found;
  int  i, max_steps;

  len = -1;

  start_mer = (char *) malloc (Kmer_Len + 1);
  strncpy (start_mer, seq + from, Kmer_Len);
  start_mer [Kmer_Len] = '\0';
  end_mer = (char *) malloc (Kmer_Len + 1);
  strncpy (end_mer, seq + to, Kmer_Len);
  end_mer [Kmer_Len] = '\0';
  fmer = (char *) malloc (Kmer_Len + 1);
  fmer [Kmer_Len] = '\0';
  rmer = (char *) malloc (Kmer_Len + 1);
  rmer [Kmer_Len] = '\0';

  max_steps = int (10 + 1.1 * (to - from));

  if (strcmp (start_mer, end_mer) == 0)
    {
      found = true;
      len = 0;
    }
  else
    {
      // Add start_mer to the used hash and to kmer_queue with a distance
      // (len field) of zero
      kmer_hash -> Kmer_To_Binary (start_mer, kref . bin);
      kref . len = 0;
      kmer_queue . push (kref);
      used . Insert (kref . bin, dummy);

      found = false;
      while (! kmer_queue . empty () && ! found)
        {
          // pop front of queue
          kref = kmer_queue . front ();
          kmer_queue . pop ();
          len = kref . len + 1;

          kmer_hash -> Binary_To_Kmer (kref . bin, fmer, true);
if (1 < Verbose)
  {
    printf ("queue size= %lu  len= %d\n", kmer_queue . size (), len);
    printf ("extracted kmer= %s\n", fmer);
  }
          strcpy (rmer, fmer);
          Reverse_Complement (rmer);
          memmove (fmer, fmer + 1, Kmer_Len - 1);  // shift left
          memmove (rmer + 1, rmer, Kmer_Len - 1);  // shift right

          // for each successor
          for (i = 0; i < 4; i ++)
            {
              fmer [Kmer_Len - 1] = FWD_ACGT [i];

              // if == endmer, set found and break; len already has the right value
              if (strcmp (fmer, end_mer) == 0)
                {
                  found = true;
                  break;
                }

              // if not in kmer_hash, skip it
              rmer [0] = REV_ACGT [i];
              p = Find_Mer (kmer_hash, fmer, rmer);
if (1 < Verbose)
  {
    printf ("fmer= %s\n", fmer);
    printf ("rmer= %s\n", rmer);
    printf ("found= %s\n", (p == NULL ? "false" : "true"));
    //**ALD left off here
  }
              if (p == NULL)
                continue;

              // if not used add to kmer_queue and to used
              // with a distance one more than front element distance
              p = Find_Mer (& used, fmer, rmer);
              if (p == NULL && len < max_steps)
                {
                  kmer_hash -> Kmer_To_Binary (fmer, kref . bin);
                  kref . len = len;
                  kmer_queue . push (kref);
                  used . Insert (kref . bin, dummy);
                }
            }
        }
    }
//**ALD Left off here
  if (found)
    printf ("Found shortest path, len= %d\n", len);
  else
    printf ("No path found\n");

  used . Clear ();

  free (start_mer);
  free (end_mer);
  free (fmer);
  free (rmer);

  return found;
}


static int  Prefix_Edit_Dist
  (const char * a, int m, const char * b, int n, int & a_hi, int & b_hi)

// Find the best edit distance match between prefixes of strings a and
// b (of lengths m & n, respectively), where the match must extend to
// the end of either a or b.  Return the edit distance and set a_hi
// and b_hi to the the length of the matched part of a and b,
// respectively.


{
  static short int  * space = NULL;
  static int  save_sz = -1;
  int  best_d, from, max_errs, mx, sz;
  int  d, e, i, j;

  if (1 < Verbose)
    {
      printf ("Prefix_Edit_Dist:\n");
      printf ("  a= %-.20s  m= %d\n", a, m);
      printf ("  b= %-.20s  n= %d\n", b, n);
    }

  max_errs = int (3.5 + 0.15 * Min (m, n));  // Allow up to 15% error + 3
  sz = (1 + max_errs) * (1 + max_errs);

  if (save_sz < sz)
    {
      space = (short int *) SAFE_REALLOC (space, sz * sizeof (short int));
      save_sz = sz;
    }

#define  S(e,d)  space [(e) * ((e) + 1) + (d)]
  // Treat space as a pyramidal array with one cell on row 0, 3 on row 1,
  // 5 on row 2, ....  Cells on row e are indexed by d = -e ... +e

  for (i = 0; i < m && i < n && a [i] == b [i]; i ++)
    ;
  S (0, 0) = i;
  if (i == m || i == n)
    { // identical match
      a_hi = i;
      b_hi = i;
      return  0;
    }
  if (1 < Verbose)
    printf ("e= %d  d= %d  S= %d  max_errs= %d\n",
            0, 0, S (0, 0), max_errs);

  for (e = 1; e <= max_errs; e ++)
    {
      mx = 0;
      for (d = - e; d <= e; d ++)
        {
          i = -1;
          if (- e < d && d < e)
            {
              i = 1 + S (e - 1, d);
              from = d;
            }
          if (- e < d - 1 && i < (j = S (e - 1, d - 1)))
            {
              i = j;
              from = d - 1;
            }
          if (d + 1 < e && i < (j = 1 + S (e - 1, d + 1)))
            {
              i = j;
              from = d + 1;
            }
          while (i < m && i + d < n && a [i] == b [i + d])
            i ++;
          S (e, d) = i;

          if (1 < Verbose)
            printf ("e= %d  d= %d  from= %d  S= %d  sub= %d\n",
                    e, d, from, S (e, d), e * (e + 1) + d);

          if (i == m || i + d == n)
            {
              a_hi = i;
              b_hi = i + d;
              return  e;
            }

          if (mx < i)
            {
              mx = i;
              best_d = d;
            }
        }
    }

  a_hi = mx;
  b_hi = mx + best_d;

  return max_errs;
}


void  Reverse_Complement
  (char * s)

// Set string  s  to its DNA Watson-Crick reverse complement

{
  int  i, j, n;

  n = strlen (s);
  for  (i = 0, j = n - 1;  i < j;  i ++, j --)
    {
      char  ch;

      ch = s [j];
      s [j] = Complement (s [i]);
      s [i] = Complement (ch);
    }

  if  (i == j)
    s [i] = Complement (s [i]);

  return;
}


static int  Suffix_Edit_Dist
  (const char * a, int a_hi, const char * b, int & a_lo, int & b_lo)

// Find the best edit distance match between a suffix of a [0..a_hi]
// and string b, where the match must extend to the beginning of
// either a or b.  Return the edit distance and set a_lo and b_lo to
// the positions in a and b, respectively, where the match ends.

{
  static short int  * space = NULL;
  static int  save_sz = -1;
  const char  * aa, * bb;
    // point to right ends of a & b, resp, to make it easier to go left
  int  best_d, from, max_errs, mx, sz;
  int  d, e, i, j, m, n;

  if (1 < Verbose)
    printf ("Suffix_Edit_Dist:\n");
  m = a_hi + 1;
  n = strlen (b);

  aa = a + a_hi;
  bb = b + n - 1;

  max_errs = int (3.5 + 0.15 * Min (m, n));  // Allow up to 15% error + 3
  sz = (1 + max_errs) * (1 + max_errs);

  if (save_sz < sz)
    {
      space = (short int *) SAFE_REALLOC (space, sz * sizeof (short int));
      save_sz = sz;
    }

#define  S(e,d)  space [(e) * ((e) + 1) + (d)]
  // Treat space as a pyramidal array with one cell on row 0, 3 on row 1,
  // 5 on row 2, ....  Cells on row e are indexed by d = -e ... +e

  for (i = 0; i < m && i < n && aa [- i] == bb [- i]; i ++)
    ;
  S (0, 0) = i;
  if (i == m || i == n)
    { // identical match
      a_lo = m - i;
      b_lo = n - i;
      return  0;
    }
  if (1 < Verbose)
    printf ("e= %d  d= %d  S= %d\n", 0, 0, S (0, 0));

  for (e = 1; e <= max_errs; e ++)
    {
      mx = 0;
      for (d = - e; d <= e; d ++)
        {
          i = -1;
          if (- e < d && d < e)
            {
              i = 1 + S (e - 1, d);
              from = d;
            }
          if (- e < d - 1 && i < (j = S (e - 1, d - 1)))
            {
              i = j;
              from = d - 1;
            }
          if (d + 1 < e && i < (j = 1 + S (e - 1, d + 1)))
            {
              i = j;
              from = d + 1;
            }
          while (i < m && i + d < n && aa [- i] == bb [-i - d])
            i ++;
          S (e, d) = i;
          if (i == m || i + d == n)
            {
              a_lo = m - i;
              b_lo = n - i - d;
              return  e;
            }
          if (1 < Verbose)
            printf ("e= %d  d= %d  from= %d  S= %d  sub= %d\n",
                    e, d, from, S (e, d), e * (e + 1) + d);

          if (mx < i)
            {
              mx = i;
              best_d = d;
            }
        }
    }

  a_lo = m - mx;
  b_lo = n - mx - best_d;

  return max_errs;
}


static void  Usage
  (void)

// Print to stderr description of options and command line for
// this program.

{
  fprintf (stderr,
    "USAGE:  kmer-repair [options] <sequence-file> <kmer-counts>\n"
    "\n"
    "Read a list of kmers with counts, in the format produced\n"
    "by 'jellyfish dump -c' from <kmer-counts>.  Then read a DNA\n"
    "multi-fasta sequence file from <sequence-file> and replace\n"
    "regions in sequence not matched by kmers with the closest sequence\n"
    "obtained by a path of consecutive kmers connecting the kmers matching\n"
    "the flanks of the unmatched section.  Either of <sequence-file>\n"
    "or <kmer-counts> can be '-' (but not both) in which case that information\n"
    "is read from stdin.  Output goes to stdout.\n"
    "\n"
    "Options:\n"
    " -h\n"
    "    Print this message\n"
    " -V <n>\n"
    "    Set verbose level to <n>; higher values for more debugging output\n"
    " -x <n>\n"
    "    Set max number of sequences used for extension/correction to <n>\n"
    "\n");

   return;
  }


