#!/bin/awk -f
# Usage:  exon1-start-mer.awk
#   Read a fasta file and extract from each sequence the first kmer
#   starting with ATG in the main reading frame (i.e., at position ~1 mod 3).
#   Output a tag, the kmer and the letter 'R' so the file can be used
#   as input to multi-walk


BEGIN {
  Kmer_Len = 71;
    # length of kmer to output
}


/^>/ {
  Pr();
  tag = substr ($1, 2);
  next;
}

{
  s = s toupper ($1);
}


END {
  Pr();
}


function Pr ()
{
  if (s != "")
    {
      n = length (s);
      p = 0;
      for (i = 1; i < n; i += 3)
        if (substr (s, i, 3) == "ATG")
          {
            p = i;
            break;
          }
      if (0 < p && p + Kmer_Len <= n + 1)
        printf "%s  %s  R\n", tag, substr (s, p, Kmer_Len);
    }

  s = "";
}


