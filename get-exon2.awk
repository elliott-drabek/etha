#!/bin/awk -f
# Usage:  get-exon2.awk
#   Read fasta sequences from stdin and output likely exon2 subsequences from
#   them.  These are sequences that begin with Start_Pattern and end at the end of
#   the longest open reading frame.


BEGIN {
  Start_Pattern = "AGAA";
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

  i = match (s, Start_Pattern);
  if (i == 0)
    return;  # ignore this string if Start_Pattern doesn't occur
  n = length (s);
  open [0] = open [1] = open [2] = 1;
  num_open = 3;
  for (j = i; j < n - 2; j ++)
    {
      f = j % 3;
      codon = substr (s, j, 3);
      if (match (codon, /TAA|TAG|TGA/) && open [f])
        {
          open [f] = 0;
          num_open --;
          if (num_open == 0)
            break;
        }
    }

  t = substr (s, i, j + 3 - i);
  print hdr;
  n = length (t);
  for (i = 1; i <= n; i += 60)
    print substr (t, i, 60);
}


