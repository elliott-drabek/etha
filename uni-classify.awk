#!/bin/awk -f
# Usage:  uni-classify.awk  <starts> <ends>
#   Read exon start patterns from starts, and exon end patterns from <ends>
#   Then read fasta sequences from stdin (presumably from unitigs, but could
#   be anything) and for each output its tag and whether it begins with
#   a start sequences or ends with an end sequence.


BEGIN {
  if  (ARGC < 2)
    Usage_Exit();

  fp = ARGV [1];
  delete ARGV [1];

  num_starts = 0;
  while  ((getline < fp) > 0)
    start [num_starts ++] = tolower ($1);

  fp = ARGV [2];
  delete ARGV [2];

  num_ends = 0;
  while  ((getline < fp) > 0)
    end [num_ends ++] = tolower ($1);

  s = "";
}


/^>/{
  Classify();
  id = substr ($1, 2);
  next;
}


{
  s = s tolower ($1);
}


END {
  Classify();
}


function Classify ()
{
  if (s != "")
    {
      n = length (s);

      is_start = 0;
      for (i = 0; i < num_starts; i ++)
        if (match (s, start [i]) == 1)
          {
            is_start = 1;
            break;
          }

      is_end = 0;
      for (i = 0; i < num_ends; i ++)
        if (match (s, end [i]) + length (end [i]) == n + 1)
          {
            is_end = 1;
            break;
          }

      printf "%s %6d   %s    %s\n", id, n, (is_start ? "start" : "  -  "),
        (is_end ? "end" : " - ");
    }

  s = "";
}


function Usage_Exit ()
{
  print "# Usage:  uni-classify.awk  <starts> <ends>";
  print "#   Read exon start patterns from starts, and exon end patterns from <ends>";
  print "#   Then read fasta sequences from stdin (presumably from unitigs, but could";
  print "#   be anything) and for each output its tag and whether it begins with";
  print "#   a start sequences or ends with an end sequence.";

  exit;
}
