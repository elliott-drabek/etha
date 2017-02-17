#!/bin/awk -f
# Usage:  exon1-end-mer.awk
#   Read a fasta file and extract from each sequence the kmer ending
#   at the last GT but prefer an earlier AAGGT if there is one within 1,000bp.
#   Output a tag, the kmer and the letter 'L'
#   so the file can be used as input to multi-walk


BEGIN {
  Kmer_Len = 71;
    # length of kmer to output
  s = "";
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
      for (i = n - 1; Kmer_Len - 1 <= i; i --)
        if (substr (s, i, 2) == "GT")
          {
            p = i;
            break;
          }

      q = 0;
      for (i --; Kmer_Len - 4 <= i; i --)
        if (substr (s, i, 5) == "AAGGT")
          {
            q = i;
            break;
          }

      if (0 < q && p - q <= 1003)
        printf "%s  %s  L\n", tag, substr (s, q + 5 - Kmer_Len, Kmer_Len);
      else if (0 < p)
        printf "%s  %s  L\n", tag, substr (s, p + 2 - Kmer_Len, Kmer_Len);
    }

  s = "";
}


