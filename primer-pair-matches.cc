//  A. L. Delcher
//
//  File:  primer-pair-matches.cc
//
//  Last Modified:  10 Mar 2014
//
//  This program reads forward and reverse primer sequences from the
//  files named as arguments (one sequence per line, no other
//  information).  It then reads a multifasta file from stdin and
//  outputs for each sequence the regions that are bounded by exact
//  primer matches.  Output is the fasta header line for each fasta
//  sequence.  Then one line for each match (if any) consisting of:
//    lo-pos  hi-pos  len  lo-primer  hi-primer
//  where lo-pos and hi-pos are the positions (counting from 1) of the
//  matched regions (including primers); len in the length of the
//  matched region; lo-primer and hi-primer are the number and which
//  primer file (fwd or rev) the primer matching each position was.


#include  "primer-pair-matches.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static char  * Fwd_Primer_Filename = NULL;
  // Name of file with list of forward primer sequences
static char  * Rev_Primer_Filename = NULL;
  // Name of file with list of reverse primer sequences


int  main
  (int argc, char * argv [])

{
  FILE  * fp;
  vector <char *>  fwd_primer, rev_primer;
  string  seq, hdr;
  char  line [MAX_LINE];
  long int  len;
  int  num_fwd, num_rev;
  int  i, n;

  Verbose = 0;

  Parse_Command_Line (argc, argv);

  fp = File_Open (Fwd_Primer_Filename, "r");
  while  (fgets (line, MAX_LINE, fp) != NULL)
    {
      n = strlen (line);
      if (0 < n && line [n - 1] == '\n')
        line [-- n] = '\0';
      for (i = 0; i < n; i ++)
        line [i] = tolower (line [i]);
      // Maybe check for duplicates here?
      fwd_primer . push_back (strdup (line));
    }
  fclose (fp);

  num_fwd = fwd_primer . size ();

  if (0 < Verbose)
    {
      printf ("\nFwd Primers:\n");
      for (i = 0; i < num_fwd; i ++)
        printf ("%2d  %s\n", i, fwd_primer [i]);
    }

  fp = File_Open (Rev_Primer_Filename, "r");
  while  (fgets (line, MAX_LINE, fp) != NULL)
    {
      n = strlen (line);
      if (0 < n && line [n - 1] == '\n')
        line [-- n] = '\0';
      for (i = 0; i < n; i ++)
        line [i] = tolower (line [i]);
      // Maybe check for duplicates here?
      rev_primer . push_back (strdup (line));
    }
  fclose (fp);

  num_rev = rev_primer . size ();

  if (0 < Verbose)
    {
      printf ("\nRev Primers:\n");
      for (i = 0; i < num_rev; i ++)
        printf ("%2d  %s\n", i, rev_primer [i]);
    }


  while (Fasta_Read (stdin, seq, hdr))
    {
      vector <Match_t>  fwd_match, rev_match;
      Match_t  m;
      size_t  p;
      int  nf, nr;
      int  j;

      len = seq . length ();

      printf (">%s\n", hdr . c_str ());

      for (i = 0; i < len; i ++)
        seq [i] = tolower (seq [i]);

      // Find matches to the forward strand and put them in fwd_match
      for (j = 0; j < num_fwd; j ++)
        {
          for (p = 0; (p = seq . find (fwd_primer [j], p)) != string :: npos; p ++)
            {
              m . primer_dir = 'f';
              m . which = j;
              m . pos = int (p + 1);  // adjust to start-at-one coords
              fwd_match . push_back (m);
            }
        }

      for (j = 0; j < num_rev; j ++)
        {
          for (p = 0; (p = seq . find (rev_primer [j], p)) != string :: npos; p ++)
            {
              m . primer_dir = 'r';
              m . which = j;
              m . pos = int (p + 1);  // adjust to start-at-one coords
              fwd_match . push_back (m);
            }
        }

      if (0 < Verbose)
        {
          nf = fwd_match . size ();
          if (0 < nf)
            printf ("Forward matches:\n");
          for (i = 0; i < nf; i ++)
            printf ("%6d %4d %c\n", fwd_match [i] . pos, fwd_match [i] . which,
                    fwd_match [i] . primer_dir);
        }

      Reverse_Complement (seq);

      // Find matches to the reverse strand and put them in rev_match
      for (j = 0; j < num_rev; j ++)
        {
          for (p = 0; (p = seq . find (fwd_primer [j], p)) != string :: npos; p ++)
            {
              m . primer_dir = 'f';
              m . which = j;
              m . pos = int (len - p);  // position of last base on fwd strand
              rev_match . push_back (m);
            }
        }

      for (j = 0; j < num_rev; j ++)
        {
          for (p = 0; (p = seq . find (rev_primer [j], p)) != string :: npos; p ++)
            {
              m . primer_dir = 'r';
              m . which = j;
              m . pos = int (len - p);  // position of last base on fwd strand
              rev_match . push_back (m);
            }
        }

      if (0 < Verbose)
        {
          nr = rev_match . size ();
          if (0 < nr)
            printf ("Reverse matches:\n");
          for (i = 0; i < nr; i ++)
            printf ("%6d %4d %c\n", rev_match [i] . pos, rev_match [i] . which,
                    rev_match [i] . primer_dir);
        }

      // Print pairs of valid matches
      nf = fwd_match . size ();
      nr = rev_match . size ();
      for (i = 0; i < nf; i ++)
        for (j = 0; j < nr; j ++)
          if (fwd_match [i] . pos < rev_match [j] . pos)
            printf ("%6d %6d %5d  %4d %c %4d %c\n",
                    fwd_match [i] . pos, rev_match [i] . pos,
                    1 + rev_match [i] . pos - fwd_match [i] . pos,
                    fwd_match [i] . which, fwd_match [i] . primer_dir,
                    rev_match [i] . which, rev_match [i] . primer_dir);

    }

  return  0;
}


static char  Complement
  (char ch)

// Returns the DNA complement of  ch

{
  return  COMPLEMENT_TABLE [unsigned (ch)];
}


static void  Parse_Command_Line
  (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

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

  if  (errflg || optind > argc - 2)
    {
      Usage ();
      exit (EXIT_FAILURE);
    }

  Fwd_Primer_Filename = argv [optind ++];
  Rev_Primer_Filename = argv [optind ++];

  return;
}


static void  Reverse_Complement
  (string & s)

// Set string  s  to its DNA Watson-Crick reverse complement

{
  unsigned  i, j, n;

  n =  s . length ();
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


static void  Usage
  (void)

//  Print to stderr description of options and command line for
//  this program.

{
  fprintf (stderr,
    "USAGE:  primer-pair-matches [options] <fwd-primers> <rev-primers>\n"
    "\n"
    "Reads forward and reverse primer sequences from files\n"
    "<fwd-primers> and <rev-primers> (one sequence per line, no other\n"
    "information).  Then read a multifasta file from stdin and\n"
    "output for each sequence the regions that are bounded by exact\n"
    "primer matches.  Output is the fasta header line for each fasta\n"
    "sequence.  Then one line for each match (if any) consisting of:\n"
    "  lo-pos  hi-pos  len  lo-primer  hi-primer\n"
    "where lo-pos and hi-pos are the positions (counting from 1) of the\n"
    "matched regions (including primers); len in the length of the\n"
    "matched region; lo-primer and hi-primer are the number and which\n"
    "primer file (fwd or rev) the primer matching each position was.\n"
    "\n"
    "Options:\n"
    " -h\n"
    "    Print this message\n"
    " -V <n>\n"
    "    Set verbose level to <n>.  Higher values are more debugging output.\n"
    "\n");

  return;
}



