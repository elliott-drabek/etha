//  A. L. Delcher
//
//  File:  multi-walk.cc
//
//  Last Modified:  13 Jan 2014
//
//  This program reads a list of kmer counts and starting kmers
//  and finds the subgraph obtained by walking from those starting kmers
//  until the end of all open reading frames.


#include  "multi-walk.hh"


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

bool  Complete_Strings = false;
  // If set true by -C option, then output complete path strings without
  // stopping when hit a previously used kmer
char  Default_Walk_Dir = 'R';
  // Presumed direction of walk if not given explicitly
bool  Ignore_ORF = false;
  // If set true by -g option, pay no attention to open reading frames
const char  * Kmer_Count_Filename;
  // Name of file from which to read kmer counts
int  Kmer_Len = -1;
  // Length of kmers
unsigned  Max_Walks = UINT_MAX;
  // Most walks from any start kmer, set by -x option
int  Num_Steps = DEFAULT_NUM_STEPS;
  // Number of steps to take in walk
char  * Start_Param;
  // Either name of file with starts or start kmer itself depending
  // on Start_On_Command_Line
bool  Start_On_Command_Line = false;
  // If true then start from the single kmer on the command line
const char  * Stop_Pattern_Filename = NULL;
  // Name of file containing stop patterns from -s option
const char  * Used_Kmer_Filename = NULL;
  // Name of file to which to output used kmers from -u option


int  main
  (int argc, char * argv [])

{
  FILE  * kmer_count_fp, * stop_pattern_fp;
  vector <char *>  start_kmer, walk_tag, stop_pattern;
  vector <char>  walk_dir;
  std :: map <string, int>  result_mers;
  unsigned  open_frame;
  Kmer_Info_t  info;
  char  fmer [MAX_LINE], rmer [MAX_LINE];
  char  a [MAX_LINE], b [MAX_LINE], s [MAX_LINE], line [MAX_LINE];
  char  wd;
  unsigned  max_count;
  int  ct, num_open_frames;
  int  i, j, k, n;

  Verbose = 0;

  Parse_Command_Line (argc, argv);

  if (Stop_Pattern_Filename != NULL)
    {
      stop_pattern_fp = File_Open (Stop_Pattern_Filename, "r");
      while  (fgets (line, MAX_LINE, stop_pattern_fp) != NULL)
        {
          n = strlen (line);
          if (line [n - 1] == '\n')
            line [-- n] = '\0';

          for (i = 0; i < n; i ++)
            line [i] = toupper (line [i]);

          stop_pattern . push_back (strdup (line));
        }
    }

  if (Start_On_Command_Line)
    {
      Kmer_Len = strlen (Start_Param);
      for (i = 0; i < Kmer_Len; i ++)
        Start_Param [i] = toupper (Start_Param [i]);
      start_kmer . push_back (strdup (Start_Param));
      walk_dir . push_back (Default_Walk_Dir);
      walk_tag . push_back (strdup ("walk_string"));
    }
  else
    {
      FILE  * fp;

      fp = File_Open (Start_Param, "r");

      while  (fgets (line, MAX_LINE, fp) != NULL)
        {
          char  tag [MAX_LINE];
          n = sscanf (line, "%s %s %s", tag, a, b);
          if (0 < n)
            {
              if (n < 2)
                {
                  sprintf (Clean_Exit_Msg_Line,
                             "ERROR:  Bad start line--need tag and sequence\n%s",
                               line);
                      Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
                }
              start_kmer . push_back (strdup (a));
              walk_tag . push_back (strdup (tag));
              if (n == 3)
                {
                  wd = toupper (b [0]);
                  if (strchr ("LRB", wd) == NULL)
                    {
                      sprintf (Clean_Exit_Msg_Line,
                               "ERROR:  Bad direction character in %s\n%s",
                               Start_Param, line);
                      Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
                    }
                  walk_dir . push_back (b [0]);
                }
              else
                walk_dir . push_back (Default_Walk_Dir);
            }
          if (Kmer_Len < 0)
            Kmer_Len = strlen (a);
          else if (int (strlen (a)) < Kmer_Len)
            {
              sprintf (Clean_Exit_Msg_Line,
                       "ERROR:  Kmer length %d is longer than length %d of start:\n%s",
                       Kmer_Len, int (strlen (a)), line);
              Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
            }
        }
    }


  if (strcmp (Kmer_Count_Filename, "-") == 0)
    kmer_count_fp = stdin;
  else
    kmer_count_fp = File_Open (Kmer_Count_Filename, "r");

  // Create the kmer hash
  Kmer_Hash_t <Kmer_Info_t>   Kmer_Hash (Kmer_Len, 11);

  max_count = (2 << KMER_INFO_FREQ_BITS) - 1;
  printf ("max_count = %d\n", max_count);
  info . used_fwd = info . used_rc = 0;

  fprintf (stderr, "Reading kmer counts ...\n");
  i = 0;
  while (fscanf (kmer_count_fp, "%s %d", s, & ct) == 2)
    {
      
      if (max_count < unsigned (ct))
        info . freq = max_count;
      else
        info . freq = ct;

      Kmer_Hash . Insert (s, info);

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

  n = start_kmer . size ();
  for (k = 0; k < n; k ++)
    {
      char  * result_string, * result_string_space = NULL;
      vector <Stack_Entry_t>  walk_stack;
      Stack_Entry_t  this_mer, next_mer [4];
      Kmer_Info_t  * p;
      string  result_s;
      bool  backtrack, hit_stop_pattern = false, occurs_already;
      char  dir;
      const char  * s;
      unsigned  variant;
      int  result_lo, result_hi;
      int  frame, num_found, pos, start_pos, step;
      int  m, w, num_walks, extra, start_len;

      if (walk_dir [k] == 'B')
        {
          dir = 'R';
          num_walks = 2;
        }
      else
        {
          dir = walk_dir [k];
          num_walks = 1;
        }

      a [Kmer_Len] = b [Kmer_Len] = '\0';
      variant = 0;

      // Clear map of mers on result string;
      result_mers . clear ();

      for (w = 0; w < num_walks; w ++)
        {
          bool  fwd_ori;     // true iff kmer on walk is the canonical one
                             //   in the hash
          char  * last_mer;  // on the end of result_string
          
          printf ("\nWalk for  %s  dir= %c\n", walk_tag [k], dir);
          printf ("%s  Start\n", start_kmer [k]);

          walk_stack . clear ();

          start_len = strlen (start_kmer [k]);
          if (dir == 'R')
            extra = start_len - Kmer_Len;
              // in case start_kmer [k] is longer than Kmer_Len
          else
            extra = 0;
          strncpy (fmer, start_kmer [k] + extra, Kmer_Len);
          strncpy (rmer, start_kmer [k] + extra, Kmer_Len);
          fmer [Kmer_Len] = rmer [Kmer_Len] = '\0';
          Reverse_Complement (rmer);

          if (w == 0)
            {
              result_string_space
                = (char *) SAFE_REALLOC (result_string_space,
                                         2 * Num_Steps + start_len + 3);
              result_string = result_string_space + Num_Steps;
              strcpy (result_string, start_kmer [k]);
              result_lo = 0;
              result_hi = start_len;
                // Working result string is result_string [result_lo .. (result_hi - 1)]
            }

          // Use just the forward strand for open reading frames
          open_frame = 0x7;  // All 3 frames open
          pos = start_pos = start_len - 1;
            // pos is the position of the rightmost base of fmer
          if (dir == 'R')
            frame = (pos + 1) % 3;
          else
            frame = 0;
            // frame is 0, 1 or 2 determined by the position of the 1 base
            // of the codon mod 3.  Going right the codon is the last 3bp
            // of the fmer; going left, it is the first 3bp
          num_open_frames = Set_Open_Frame (open_frame, start_kmer [k], pos);

          if (0 < Verbose)
            printf ("num_open_frames= %d  %c %c %c\n", num_open_frames,
                    ((open_frame & 0x1) ? 'T' : 'F'),
                    ((open_frame & 0x2) ? 'T' : 'F'),
                    ((open_frame & 0x4) ? 'T' : 'F'));

          if (strcmp (fmer, rmer) < 0)
            s = fmer;
          else
            s = rmer;

          p = Kmer_Hash . Find (s);
          if (p == NULL)
            {
              printf ("  Start kmer not found\n");
              num_found = 0;
            }
          else
            {
              printf ("Start kmer freq= %d\n", p -> freq);
              num_found = 1;
              this_mer . pos = pos;
              this_mer . r_lo = result_lo;
              this_mer . r_hi = result_hi;
              this_mer . freq = p -> freq;
              this_mer . dir = dir;
              this_mer . fmer = fmer;
              this_mer . rmer = rmer;
              this_mer . frame = frame;
              this_mer . open_frame = open_frame;
              walk_stack . push_back (this_mer);
            }
          
          backtrack = false;

          while (! walk_stack . empty () && variant < Max_Walks)
//          for (i = 0; i < Num_Steps && 0 < num_open_frames && 0 < num_found; i ++)
            {
              string  ms;
              int  pattern_sub;

              this_mer = walk_stack . back ();
              walk_stack . pop_back ();

              if (backtrack)
                { // remove eliminated kmers from result_mers
                  int  j;

                  if (dir == 'R')
                    {
                      for (j = pos; this_mer . pos <= j; j --)
                        {
                          strncpy (a, result_string + j + 1 - Kmer_Len, Kmer_Len);
                          ms = a;
                          result_mers . erase (ms);
                        }
                    }
                  else
                    {
                      for (j = pos; j <= this_mer . pos; j ++)
                        {
                          strncpy (a, result_string + j + 1 - Kmer_Len, Kmer_Len);
                          ms = a;
                          result_mers . erase (ms);
                        }
                    }
                }
              pos = this_mer . pos;
              result_lo = this_mer . r_lo;
              result_hi = this_mer . r_hi;
              result_string [result_hi] = '\0';
              frame = this_mer . frame;
              open_frame = this_mer . open_frame;
              num_open_frames = Num_Open (open_frame);
              step = abs (pos - start_pos);

              strcpy (fmer, this_mer . fmer . c_str ());
              strcpy (rmer, this_mer . rmer . c_str ());

              // Put correct end character on result string
              // result_string is built from fmer's; rmer's are for hash lookups
              if (dir == 'R')
                {
                  result_string [result_hi - 1] = fmer [Kmer_Len - 1];
                  last_mer = fmer;
                }
              else
                {
                  result_string [result_lo] = fmer [0];
                  last_mer = rmer;
                }

              if (strcmp (fmer, rmer) < 0)
                {
                  s = this_mer . fmer . c_str ();
                  fwd_ori = true;
                }
              else
                {
                  s = this_mer . rmer . c_str ();
                  fwd_ori = false;
                }
              p = Kmer_Hash . Find (s);
              assert (p != NULL);

              if (0 < Verbose)
                {
                  if (backtrack)
                    printf ("\nBacktrack to fmer= %s\n", fmer);
                  printf ("\nstep= %d  pos= %d  num_open_frames= %d  frame= %d\n"
                          "  used_fwd/rc= %c/%c  backtrack= %c\n",
                          step, pos, num_open_frames, frame,
                          Bool_To_Char (p -> used_fwd),
                          Bool_To_Char (p -> used_rc), Bool_To_Char (backtrack));
                  printf ("result_lo/hi= %d/%d  hi_ch= %c\n", result_lo,
                          result_hi, result_string [result_hi - 1]);
                  if (1 < Verbose)
                    printf ("fmer= %s  rmer= %s\n", fmer, rmer);
                }

              // Check already used first.  If so and backtracking then
              // don't need to print anything.  But if not backtracking
              // need to print

              if (! Complete_Strings
                  && ((fwd_ori && p -> used_fwd)
                      || (! fwd_ori && p -> used_rc)))
                {
                  if (0 < Verbose)
                    {
                      printf ("  already used:  ");
                      if (backtrack)
                        printf ("backtracking--no need to print\n");
                      else
                        {
                          printf ("not backtracking--print partial string\n");
                          printf ("result_lo= %d  result_hi= %d\n", result_lo,
                                  result_hi);
                        }
                    }
                  if (! backtrack)
                    Print_Result (stdout, result_string, result_lo, result_hi,
                                  walk_tag [k], ++ variant, "Hit-used-kmer");

                  backtrack = true;
                  continue;
                }

              // See if at last step, or no open reading frames, or hit
              // stop pattern.
              // If so, do nothing and loop will pop next stack value;
              // otherwise, continue.  Set a flag for which we're doing.

              if (Stop_Pattern_Filename != NULL)
                hit_stop_pattern = Check_Stops (stop_pattern, fmer, dir, pattern_sub);

              // Check for cycles by seeing if last kmer occurs earlier
              // on the path.  If so print and backtrack.
              occurs_already
                = (result_mers . find (this_mer . fmer) != result_mers . end ());
              
              if (step == Num_Steps || (num_open_frames == 0 && ! Ignore_ORF)
                    || hit_stop_pattern || occurs_already)
                {
                  char  expl [MAX_LINE];

                  if (hit_stop_pattern)
                    sprintf (expl, "Hit-stop-pattern %s freq=%d",
                             stop_pattern [pattern_sub], p -> freq);
                  else if (num_open_frames == 0 && ! Ignore_ORF)
                    sprintf (expl, "Orf-end freq=%d", p -> freq);
                  else if (step == Num_Steps)
                    sprintf (expl, "Max-steps freq=%d", p -> freq);
                  else if (occurs_already)
                    sprintf (expl, "Cycle-detected  pos= %d",
                             result_mers [this_mer . fmer]);
                  else
                    {
                      fprintf (stderr, "ERROR:  can't get here\n");
                      exit (EXIT_FAILURE);
                    }
                      
                  printf ("expl= %s\n", expl);
                  printf ("result_lo= %d  result_hi= %d\n", result_lo,
                          result_hi);
                  Print_Result (stdout, result_string, result_lo, result_hi,
                                walk_tag [k], ++ variant, expl);

                  if (fwd_ori)
                    p -> used_fwd = true;
                  else
                    p -> used_rc = true;
                  backtrack = true;
                  continue;
                }

              // Add fmer to result_mers
              result_mers [this_mer . fmer] = pos;

              if (fwd_ori)
                p -> used_fwd = true;
              else
                p -> used_rc = true;
              if (dir == 'R')
                {
                  strncpy (a, fmer + 1, Kmer_Len - 1);
                  strncpy (b + 1, rmer, Kmer_Len - 1);
                  pos ++;
                  result_string [++ result_hi] = '\0';
                  frame = (frame + 1) % 3;
                }
              else
                {
                  strncpy (a + 1, fmer, Kmer_Len - 1);
                  strncpy (b, rmer + 1, Kmer_Len - 1);
                  pos --;
                  result_lo --;
                  frame = (frame + 2) % 3;  // + 2 is the same as - 1 mod 3
                }

              // consider all possibilites for the next letter
              num_found = 0;
              for (j = 0; j < 4; j ++)
                {
                  unsigned  o_f;  // local copy of open_fram

                  o_f = open_frame;
                  if (dir == 'R')
                    {
                      a [Kmer_Len - 1] = FWD_ACGT [j];
                      b [0] = REV_ACGT [j];
                      if (Is_Stop (a + Kmer_Len - 3))
                        o_f &= ~(0x1 << frame);
                    }
                  else
                    {
                      a [0] = FWD_ACGT [j];
                      b [Kmer_Len - 1] = REV_ACGT [j];
                      if (Is_Stop (a))
                        o_f &= ~(0x1 << frame);
                    }
                  if (strcmp (a, b) < 0)
                    {
                      s = a;
                      fwd_ori = true;
                    }
                  else
                    {
                      s = b;
                      fwd_ori = false;
                    }
//                  printf ("\nj= %d  a= %s  b= %s\n", j, a, b);

                  p = Kmer_Hash . Find (s);
                  if (p == NULL)
                    ;  // ignore
                  else
                    {
//                      printf ("\nFound %s  ct= %d  pos= %d\n", s, p -> freq, pos);
                      printf ("\nFound %s  ct= %d  pos= %d\n", a, p -> freq, pos);
                      if ((fwd_ori && p -> used_fwd)
                          || (! fwd_ori && p -> used_rc))
                        printf ("  already used\n");

                      // Include it even if it's already used.  Will detect
                      // that when pop from the stack
                      next_mer [num_found] . pos = pos;
                      next_mer [num_found] . r_lo = result_lo;
                      next_mer [num_found] . r_hi = result_hi;
                      next_mer [num_found] . freq = p -> freq;
                      next_mer [num_found] . dir = dir;
                      next_mer [num_found] . fmer = a;
                      next_mer [num_found] . rmer = b;
                      next_mer [num_found] . frame = frame;
                      next_mer [num_found] . open_frame = o_f;
                      num_found ++;
                    }
                }

              printf ("num_found= %d", num_found);
              if (num_found < 2)
                  printf ("\n");
              else
                  printf ("  branch_pos= %d\n", result_hi - result_lo);
 

              if (num_found == 0)
                {
                  printf ("No extension found--print partial orf\n"
                          "  #open-frames= %d  pos= %d\n", num_open_frames, pos);
                  // Set result_lo/hi back to previous value for printing
                  if (dir ==  'R')
                    result_hi --;
                  else
                    result_lo ++;
                  Print_Result (stdout, result_string, result_lo, result_hi,
                                walk_tag [k], ++ variant, "No-extension");
                  backtrack = true;
                }
              else
                {
                  // sort entries by freq
                  for (j = 0; j < num_found - 1; j ++)
                    for (m = j + 1; m < num_found; m ++)
                      if (next_mer [m] . freq < next_mer [j] . freq)
                        {
                          Stack_Entry_t  tmp = next_mer [j];

                          next_mer [j] = next_mer [m];
                          next_mer [m] = tmp;
                        }

                  // Push all possibilities onto stack
                  for (j = 0; j < num_found; j ++)
                    {
                      if (1 < Verbose)
                        printf ("push fmer= %s  freq= %d\n",
                                next_mer [j] . fmer . c_str (),
                                next_mer [j] . freq);
                      walk_stack . push_back (next_mer [j]);
                    }
                  backtrack = false;
                }
            }
          if (0 < Verbose)
            printf ("\nWalk stack has %d entries\n", int (walk_stack . size ()));

          dir = 'L';
        }

#if 0
      // print resulting string
      printf ("\n");
      for (j = result_lo; j < result_hi; j ++)
        result_string [j] = toupper (result_string [j]);
      Fasta_Print (stdout, result_string, walk_tag [k]);
      printf ("\n");
#endif

    }

  if (kmer_count_fp != stdin)
    fclose (kmer_count_fp);

  if (Used_Kmer_Filename != NULL)
    {
      FILE  * fp = File_Open (Used_Kmer_Filename, "w");

      Kmer_Hash . Dump_Kmers_Select (fp, & Print_If_Used);
    }

  return  0;
}


static bool  Check_Stops
  (vector <char *> & stop_pattern, char * fmer, char dir, int & pattern_sub)

// Check whether fmer contains any pattern in stop_pattern on the end
// corresponding to dir.  If so, return true and set pattern_sub to the
// subscript of the pattern; otherwise, return false and set pattern_sub
// to -1.

{
  int  i, m, n, offset;

  n = stop_pattern . size ();

  for (i = 0; i < n; i ++)
    {
      m = strlen (stop_pattern [i]);
      if (dir == 'R')
        offset = Kmer_Len - m;
      else
        offset = 0;
      if (strncmp (fmer + offset, stop_pattern [i], m) == 0)
        {
          pattern_sub = i;
          return true;
        }
    }

  pattern_sub = -1;
  return false;
}


char  Complement
  (char ch)

// Returns the DNA complement of  ch

{
  return  COMPLEMENT_TABLE [unsigned (ch)];
}


static bool  Is_Stop
  (const char * codon)

// Return true iff the first 3 bases of codon are a stop codon
// Assume all characters are uppercase

{
  return (strncmp (codon, "TAA", 3) == 0
          || strncmp (codon, "TAG", 3) == 0
          || strncmp (codon, "TGA", 3) == 0);
}


static int  Num_Open
  (unsigned char open_frame)

// Return the number of open reading frames in bit pattern open_fram

{
  unsigned char  mask = 0x1;
  int  i, n;

  n = (mask & open_frame);

  for (i = 0; i < 2; i ++)
    {
      mask <<= 1;
      n += ((mask & open_frame) ? 1 : 0);
    }

  return  n;
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
          && ((ch = getopt (argc, argv, "bCghik:n:rs:u:V:x:")) != EOF))
    switch  (ch)
      {
      case  'b' :
        Default_Walk_Dir = 'B';
        break;

      case  'C' :
        Complete_Strings = true;
        break;

      case  'g' :
        Ignore_ORF = true;
        break;

      case  'h' :
        Usage ();
        exit (EXIT_SUCCESS);

      case  'i' :
        Start_On_Command_Line = true;
        break;

      case  'k' :
        Kmer_Len = strtol (optarg, NULL, 10);
        break;

      case  'n' :
        Num_Steps = strtol (optarg, NULL, 10);
        break;

      case  'r' :
        Default_Walk_Dir = 'R';
        break;

      case  's' :
        Stop_Pattern_Filename = optarg;
        break;

      case  'u' :
        Used_Kmer_Filename = optarg;
        break;

      case  'V' :
        Verbose = strtol (optarg, NULL, 10);
        break;

      case  'x' :
        Max_Walks = strtol (optarg, NULL, 10);
        break;

      default :
        errflg = true;
      }

  if  (errflg || optind > argc - 2)
    {
      Usage ();
      exit (EXIT_FAILURE);
    }

  Start_Param = argv [optind ++];
  Kmer_Count_Filename = argv [optind ++];

  return;
}


void  Print_If_Used
  (FILE * fp, const char * s, const Kmer_Info_t & info)

// Print string s and info . freq to fp if either info . used_fwd
// or info . used_rc is true.

{
  if (info . used_fwd || info . used_rc)
    fprintf (fp, "%s %6d\n", s, info . freq);

  return;
}


static void  Print_Result
  (FILE * fp, char * s, int lo, int hi, char * id, unsigned ver, const char * expl)

// Print to fp the contents of s from positions lo .. (hi - 1) in fasta
// format.  Use id and the version number ver as the fasta-id and
// also include expl on the fasta header line.  Note that this function
// sets s [hi] to '\0';

{
  char hdr [MAX_LINE];

  sprintf (hdr, "%s-%03u  %s", id, ver, expl);
  s [hi] = '\0';
  fprintf (fp, "\n");
  Fasta_Print (fp, s + lo, hdr);
  fprintf (fp, "\n");

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


static int  Set_Open_Frame
  (unsigned char open_frame, const char * s, int offset)

// Set values in open_frame false iff there is a stop codon in
// the corresponding reading frame in s.  offset is the position
// of the *LAST* base in s (in a presumed longer string) to determine
// the reference open reading frame and offset may be negative.  Return
// the number of open reading frames (i.e., trues) in open_frame.
// In open_frame, use 1 bit for each frame: 0x1 for frame 0; 0x2 for frame 1;
// 0x4 for frame 2.  Or generally, (1 << f) for frame f

{
  char  codon [4];
  int  f, i, k, n;

  // Reading frame 0 is when first base of codon is at position ..., -3, 0, 3, ...
  // Reading frame 1 is when first base of codon is at position ..., -2, 1, 4, ...
  // Reading frame 2 is when first base of codon is at position ..., -1, 2, 5, ...

  n = strlen (s);
  f = (1 + offset - n);  // ref position of 1st base of s
  if (f < 0)
    f += 3 * ((3 - f) % 3);  // get positive equivalent
  f %= 3;

  codon [3] = '\0';
  for (i = 0; i < n - 2; i ++)
    {
      strncpy (codon, s, 3);
      if (Is_Stop (codon))
        open_frame &= (~(0x1 << f));
      f = (f + 1) % 3;
    }

  for (i = k = 0; i < 3; i ++)
    if (open_frame & (0x1 << i))
      k ++;

  return k;
}


static void  Usage
  (void)

// Print to stderr description of options and command line for
// this program.

{
  fprintf (stderr,
    "USAGE:  multi-walk [options] <start-kmers> <kmer-counts>\n"
    "\n"
    "Read a list of kmer counts (e.g., as produced by jellyfish) from\n"
    "file <kmer-counts> and list of start kmers from file <start-kmers>.\n"
    "Find and output the subgraph obtained by walking consecutive kmers\n"
    "from the starting kmers up until a stop codon is hit in all reading\n"
    "frames.  Format of <start-kmers> is one entry per line:\n"
    "  <tag>  <kmer>  <walk direction>\n"
    "where walk direction is optional and one of:  'R' right (5' to 3'),\n"
    "'L' left (3' to 5'), or 'B' both directions\n"
    "\n"
    "Options:\n"
    " -b\n"
    "    Walk in both directions.  Applies if -i option or no direction\n"
    "    is given in <start-kmers> file\n"
    " -C\n"
    "    Output complete path strings without stopping when a previously\n"
    "    used kmer is hit\n"
    " -g\n"
    "    Ignore stop codons and walk through them\n"
    " -h\n"
    "    Print this message\n"
    " -i\n"
    "    <start-kmers> is the actual start kmer instead of a file\n"
    " -k <n>\n"
    "    Use kmer length <n> if <start-mers> has sequences longer than <n>\n"
    " -n <n>\n"
    "    Walk for at most <n> steps\n"
    " -r\n"
    "    Walk in reverse direction.  Applies if -i option or no direction\n"
    "    is given in <start-kmers> file\n"
    " -s <f>\n"
    "    Read patterns from file <f> (one per line) and stop walk when\n"
    "    any pattern is hit\n"
    " -u <f>\n"
    "    Output all kmers used in walks to file <f>\n"
    " -x <n>\n"
    "    Allow at most <n> walks from any start kmer\n"
    "\n");

   return;
  }


