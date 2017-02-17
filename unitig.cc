//  A. L. Delcher
//
//  File:  unitig.cc
//
//  Last Modified:  14 Oct 2013
//
//  This program reads (from stdin) output from
//  'show-coords -cHlT | sort -k8n -k12,12 -k1n -k2nr' from a NUCmer
//  of a set of sequences against themselves and then
//  computes unitigs for the sequences.  Output is the layout of
//  sequences in each unitig.  Only matches in the positive orientation
//  are used.


#include  "unitig.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

const int  Dummy_Constant = 0;
  // Place holder


int  main
  (int argc, char * argv [])

{
  vector <Seq_t>  seq;
  Olap_t  new_olap;
  int  a_lo, a_hi, b_lo, b_hi, dummy, a_len, b_len, a_sub, b_sub;
  double  idnty, a_frac, b_frac;
  char  a_id [MAX_LINE], b_id [MAX_LINE], prev_id [MAX_LINE];
  int  c, i, n;

  Verbose = 0;

  Parse_Command_Line (argc, argv);

  prev_id [0] = '\0';  // empty string

  // Read input
  n = 0;
  while  (scanf ("%d %d %d %d %d %d %lf %d %d %lf %lf %s %s",
                 & a_lo, & a_hi, & b_lo, & b_hi, & dummy, & dummy,
                 & idnty, & a_len, & b_len, & a_frac, & b_frac,
                 a_id, b_id) == 13)
    {
      if (strcmp (a_id, prev_id) != 0)
        {
          a_sub = Find_Or_Insert (a_id, a_len, seq);
          strcpy (prev_id, a_id);
        }

      if (b_hi < b_lo)
        continue;  // ignore negative orientations
      if (a_lo == 1 && a_hi == a_len
            && (a_len < b_len || strcmp (a_id, b_id) < 0))
               // if exact duplicate keep higher id
        {
          seq [a_sub] . contained = 1;  // mark as contained
          continue;  // skip contained sequences
        }
      if (strcmp (a_id, b_id) == 0)
        continue;  // skip self matches
      if ((a_lo != 1 && b_lo != 1) || (a_hi != a_len && b_hi != b_len))
        continue;  // skip non-overlap matches
      n ++;

      if (1 < Verbose)
        printf ("a_id= %s  b_id= %s  prev_id= %s\n", a_id, b_id, prev_id);

      b_sub = Find_Or_Insert (b_id, b_len, seq);

      new_olap . a_lo = a_lo;
      new_olap . a_hi = a_hi;
      new_olap . b_lo = b_lo;
      new_olap . b_hi = b_hi;
      new_olap . a_sub = a_sub;
      new_olap . b_sub = b_sub;

      if (a_lo == 1)
        seq [a_sub] . left_olap . push_back (new_olap);
      else
        seq [a_sub] . right_olap . push_back (new_olap);

      if (1 < Verbose)
        printf ("a_sub= %d  b_sub= %d\n", a_sub, b_sub);
    }

  printf ("%d lines used\n", n);
  n = seq . size ();
  printf ("%d different sequences\n", n);
  for (c = i = 0; i < n; i ++)
    if (seq [i] . contained)
      c ++;
  printf ("  %d are contained\n", c);

  Remove_Contained_Olaps (seq);

  Remove_Transitive_Olaps (seq);

  //**ALD left off here

  Find_Unitigs (seq);
  

  return  0;
}


static int  Find_Or_Insert
  (const char * b_id, int b_len, vector <Seq_t> & seq)

// Search for b_id in seq.  If found return its subscript; otherwise
// append an entry for it (including its length b_len) to seq and return
// that subscript.

{
  Seq_t  s;
  int  i, n;

  n = seq . size ();
  for (i = 0; i < n; i ++)
    if (strcmp (b_id, seq [i] . id) == 0)
      return i;

  s . id = strdup (b_id);
  s . len = b_len;
  seq . push_back (s);

  return n;
}

  
static void  Find_Unitigs
  (vector <Seq_t> & seq)

// Find and output unitigs comprised of non-contained sequences in seq.

{
  int  unitig_num = 0;
  int  i, n;

  n = seq . size ();
  for (i = 0; i < n; i ++)
    if (! seq [i] . contained && ! seq [i] . used)
      {
        bool  go_left, go_right;
        int  j, k;

        // Follow mutually unique overlaps to left boundary
        j = i;
        go_left = true;
        do
          {
            if (seq [j] . left_olap . size () != 1)
              {
                go_left = false;
                break;
              }
            k = seq [j] . left_olap [0] . b_sub;
            if (seq [k] . contained || seq [k] . used
                || seq [k] . right_olap . size () != 1)
              {
                go_left = false;
                break;
              }
            j = k;
          } while (go_left);

        printf ("\nUnitig #%d\n", ++ unitig_num);
        printf ("  %3d %5d %5d   %5d  %s\n",
                j, 1, seq [j] . len, 1, seq [j] . id);
        seq [j] . used = 1;
        go_right = true;
        do
          {
            if (seq [j] . right_olap . size () != 1)
              {
                go_right = false;
                break;
              }
            k = seq [j] . right_olap [0] . b_sub;
            if (seq [k] . contained || seq [k] . used
                || seq [k] . left_olap . size () != 1)
              {
                go_right = false;
                break;
              }
            printf ("  %3d %5d %5d   %5d  %s\n",
                    k, seq [j] . right_olap [0] . b_hi + 1,
                    seq [k] . len, seq [j] . right_olap [0] . a_lo,
                    seq [k] . id);
            j = k;
            seq [j] . used = 1;
          } while (go_right);
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
          && ((ch = getopt (argc, argv, "hV:")) != EOF))
    switch  (ch)
      {
      case  'h' :
        Usage ();
        exit (EXIT_SUCCESS);

      case  'V' :
        Verbose = strtol (optarg, NULL, 10);
        break;

      default :
        errflg = true;
      }

  if  (errflg || optind > argc)
    {
      Usage ();
      exit (EXIT_FAILURE);
    }

  return;
}


static void  Remove_Contained_Olaps
  (vector <Seq_t> & seq)

// Remove olaps from non-contained sequences to contained sequences in seq

{
  int  i, j, k, m, n;

  n = seq . size ();
  for (i = 0; i < n; i ++)
    if (! seq [i] . contained)
      {  // ignore olaps of contained sequences
        m = seq [i] . left_olap . size ();
        for (j = k = 0; j < m; j ++)
          {
            int  p = seq [i] . left_olap [j] . b_sub;

            if (! seq [p] . contained)
              {
                if (k < j)
                  seq [i] . left_olap [k] = seq [i] . left_olap [j];
                k ++;
              }
          }
        if (k < m)
          seq [i] . left_olap . resize (k);

        m = seq [i] . right_olap . size ();
        for (j = k = 0; j < m; j ++)
          {
            int  p = seq [i] . right_olap [j] . b_sub;

            if (! seq [p] . contained)
              {
                if (k < j)
                  seq [i] . right_olap [k] = seq [i] . right_olap [j];
                k ++;
              }
          }
        if (k < m)
          seq [i] . right_olap . resize (k);
      }
}


static void  Remove_Transitive_Olaps
  (vector <Seq_t> & seq)

// Remove from seq overlaps of the form AC if overlaps AB and BC are
// present

{
  int  i, j, k, m, n, p, q;

  n = seq . size ();

  // overlaps should be in order by amount of overlap; largest to smallest
  for (i = 0; i < n; i ++)
    {
      // skip contained sequences
      if (seq [i] . contained)
        continue;

      // do right overlaps first
      m = seq [i] . right_olap . size ();
      if (0 < Verbose)
        printf ("seq %d has %d right olaps\n", i, m);

      if (1 < m)
        {
          k = seq [i] . right_olap [0] . b_sub;
          for (j = q = 1; j < m; j ++)
            {
              p = seq [i] . right_olap [j] . b_sub;
              if (0 < Verbose)
                printf ("  overlap j= %d to %3d  len= %4d  ", j, p,
                        1 + seq [i] . len - seq [i] . right_olap [j] . a_lo);

              // have overlap i to k and i to p
              // if there is an overlap from k to p, eliminate i to p overlap
              if (seq [k] . Has_Right_Olap_To (p))
                {
                  if (0 < Verbose)
                    printf ("transitively removed\n");
                }
              else
                {
                  if (0 < Verbose)
                    printf ("kept\n");
                  if (q < j)
                    seq [i] . right_olap [q] = seq [i] . right_olap [j];
                  q ++;
                }
            }
          if (q < m)
            {
              seq [i] . right_olap . resize (q);
              if (0 < Verbose)
                printf ("  resized to %d overlaps\n", q);
            }
        }

      // now do left overlaps
      m = seq [i] . left_olap . size ();
      if (0 < Verbose)
        printf ("seq %d has %d left olaps\n", i, m);

      if (1 < m)
        {
          k = seq [i] . left_olap [0] . b_sub;
          for (j = q = 1; j < m; j ++)
            {
              p = seq [i] . left_olap [j] . b_sub;
              if (0 < Verbose)
                printf ("  overlap j= %d to %3d  len= %4d  ", j, p,
                        seq [i] . left_olap [j] . a_hi);

              // have overlap i to k and i to p
              // if there is an overlap from k to p, eliminate i to p overlap
              if (seq [k] . Has_Left_Olap_To (p))
                {
                  if (0 < Verbose)
                    printf ("transitively removed\n");
                }
              else
                {
                  if (0 < Verbose)
                    printf ("kept\n");
                  if (q < j)
                    seq [i] . left_olap [q] = seq [i] . left_olap [j];
                  q ++;
                }
            }
          if (q < m)
            {
              seq [i] . left_olap . resize (q);
              if (0 < Verbose)
                printf ("  resized to %d overlaps\n", q);
            }
        }
    }

  return;
}


static void  Usage
  (void)

// Print to stderr description of options and command line for
// this program.

{
  fprintf (stderr,
    "USAGE:  unitig < show-coords-output"
    "\n"
    "This program reads (from stdin) output from\n"
    "'show-coords -cHlT | sort -k12,12 -k1n -k2nr' from a NUCmer\n"
    "of a set of sequences against themselves and then\n"
    "computes unitigs for the sequences.  Output is the layout of\n"
    "sequences in each unitig.  Only matches in the positive orientation\n"
    "are used.\n"
    "\n"
    "Options:\n"
    " -h\n"
    "    Print this message\n"
    " -V <n>\n"
    "    Set verbose level to <n>; higher values for more debugging output\n"
    "\n");

   return;
  }


