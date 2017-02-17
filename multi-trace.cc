//  A. L. Delcher
//
//  File:  multi-trace.cc
//
//  Last Modified:  3 May 2013
//
//  This program reads a list of kmers with counts, in the format produced
//  by 'jellyfish dump -c'.  It then reads a DNA multi-fasta sequence file
//  and prints the count of each kmer in the sequences.  It also checks for
//  and reports any (k-1)mers in the sequence with a different preceding or
//  trailing base than is in the sequence that have a count entry.


#include  "multi-trace.hh"


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
const char  * Sequence_Filename;
  // Name of file with DNA fasta sequences


int  main
  (int argc, char * argv [])

{
  FILE  * sequence_fp, * kmer_count_fp;
  Kmer_Hash_t <Kmer_Info_t>  * kmer_hash;
  string  seq_string, seq_hdr;
  Kmer_Info_t  info;
  char  fmer [MAX_LINE], rmer [MAX_LINE];
  char  a [MAX_LINE], b [MAX_LINE], s [MAX_LINE];
  const char  * seq;
  unsigned  max_count;
  int  ct, seq_len;
  int  i, j;

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

      printf ("\n>%s  len=%d\n", seq_hdr . c_str (), seq_len);

      strncpy (fmer, seq, Kmer_Len);
      fmer [Kmer_Len] = '\0';
      strcpy (rmer, fmer);
      Reverse_Complement (rmer);
      
      for (i = 1; i <= seq_len - Kmer_Len + 1; i ++)
        {
          Kmer_Info_t  * p;

          p = Find_Mer (kmer_hash, fmer, rmer);
          if (p == NULL)
            ct = 0;
          else
            {
              ct = p -> freq;
              p -> hit = 1;
            }

          printf ("%5d  %s  %3d", i, fmer, ct);

          // Look for alternate left characters
          strcpy (a, fmer);
          strcpy (b, rmer);
          for (j = 0; j < 4; j ++)
            if (fmer [0] != FWD_ACGT [j])
              {
                a [0] = FWD_ACGT [j];
                b [Kmer_Len - 1] = REV_ACGT [j];
                p = Find_Mer (kmer_hash, a, b);
                if (p != NULL)
                  printf ("  %c,%c,%-5d", 'L', a [0], p -> freq);
              }

          // Look for alternate right characters
          a [0] = fmer [0];
          b [Kmer_Len - 1] = rmer [Kmer_Len - 1];
          for (j = 0; j < 4; j ++)
            if (fmer [Kmer_Len - 1] != FWD_ACGT [j])
              {
                a [Kmer_Len - 1] = FWD_ACGT [j];
                b [0] = REV_ACGT [j];
                p = Find_Mer (kmer_hash, a, b);
                if (p != NULL)
                  printf ("  %c,%c,%-5d", 'R', a [Kmer_Len - 1], p -> freq);
              }
          printf ("\n");

          // Advance to next character
          a [Kmer_Len - 1] = fmer [Kmer_Len - 1];
          b [0] = rmer [0];
          strncpy (fmer, a + 1, Kmer_Len - 1);
          strncpy (rmer + 1, b, Kmer_Len - 1);
          fmer [Kmer_Len - 1] = seq [i + Kmer_Len - 1];
          rmer [0] = Complement (fmer [Kmer_Len - 1]);
        }
    }

  if (sequence_fp != stdin)
    fclose (sequence_fp);

  fprintf (stderr, "Total kmers in hash:  %d\n", kmer_hash -> Count_All ());
  fprintf (stderr, " which were matched:  %d\n",
           kmer_hash -> Count_All_Select (& Has_Positive_Hit));

  return  0;
}


char  Complement
  (char ch)

// Returns the DNA complement of  ch

{
  return  COMPLEMENT_TABLE [unsigned (ch)];
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


bool  Has_Positive_Hit
  (const Kmer_Info_t & info)

// Return true iff info . hit > 0

{
  return (0 < info . hit);
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
          && ((ch = getopt (argc, argv, "bhik:n:rV:")) != EOF))
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

  Sequence_Filename = argv [optind ++];
  Kmer_Count_Filename = argv [optind ++];

  return;
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


static void  Usage
  (void)

// Print to stderr description of options and command line for
// this program.

{
  fprintf (stderr,
    "USAGE:  multi-trace [options] <sequence-file> <kmer-counts>\n"
    "\n"
    "Read a list of kmers with counts, in the format produced\n"
    "by 'jellyfish dump -c' from <kmer-counts>.  Then read a DNA\n"
    "multi-fasta sequence file from <sequence-file> and print the\n"
    "count of each kmer in the sequences.  Also check for and report\n"
    "any (k-1)mers in the sequence with an alternative preceding or\n"
    "trailing base from what is in the sequence and that have a count entry.\n"
    "Sequences are assumed to be non-circular.  Either of <sequence-file>\n"
    "or <kmer-counts> can be '-' (but not both) in which case that information\n"
    "is read from stdin.  Output goes to stdout.\n"
    "\n"
    "Options:\n"
    " -h\n"
    "    Print this message\n"
    " -V <n>\n"
    "    Set verbose level to <n>; higher values for more debugging output\n"
    "\n");

   return;
  }


