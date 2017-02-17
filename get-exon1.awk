#!/bin/awk -f
# Usage:  get-exon1.awk
#   Read fasta sequences from stdin and output likely exon1 subsequences from
#   them.  These are sequences that have End_Pattern near the end and from
#   that point extend back to the end of all reading frames.  In the longest
#   open reading frame, the first ATG is taken as the start.  If there is no
#   End_Pattern near the end of the sequence, the extracted sequence extends
#   to the end.  Similarly, if there is an open reading frame that extends
#   to the front of the sequence, the subsequence begins at the start of the
#   sequence.


BEGIN {
  End_Pattern = "AAGGT";
  End_Region_Len = 10;
    # length of region at end in which to look for End_Pattern
  s = "";
}


/^>/ {
  Pr();
  hdr = $0;
  s = "";
  next;
}


{
  s = s toupper ($0);
}


END {
  Pr();
}


function Pr ()
{
  if (s == "")
    return;

  last = n = length (s);
  end_pat_len = length (End_Pattern);


  for (i = n + 1 - end_pat_len; n - End_Region_Len < i; i --)
    if (substr (s, i, end_pat_len) == End_Pattern)
      {
        last = i + end_pat_len - 1;
        break;
      }

  open [0] = open [1] = open [2] = 1;
  num_open = 3;
  for (i = last - 2; 0 < i; i --)
    {
      f = i % 3;
      codon = substr (s, i, 3);
      if (match (codon, /TAA|TAG|TGA/) && open [f])
        {
          open [f] = 0;
          num_open --;
          if (num_open == 0)
            break;
        }
    }

  if (0 < num_open)
    first = 1;
  else
    {
      for (i = 1; i < last - 2; i += 3)
        if (substr (s, i, 3) == "ATG")
          {
            first = i;
            break;
          }
    }

  if (last < first)
    return;  # nothing to print

  t = substr (s, first, 1 + last - first);
  print hdr;
  n = length (t);
  for (i = 1; i <= n; i += 60)
    print substr (t, i, 60);
}


