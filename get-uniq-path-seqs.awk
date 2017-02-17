#!/bin/awk -f
# Usage:  get-uniq-path-seqs.awk
#   Read ex1.cand.fa file and output sequences that represent the only
#   extension possible from the start kmer.  These are sequences whose
#   fasta tag ends in -001 and for which there is no -002 sequence.


BEGIN {
  n = 0;
  prev = "";
}


/^>/ {
  a = match ($1, /\-[0-9]*$/);
  suffix = substr ($1, a);
  if (suffix == "-001" && prev == "-001")
    Pr();
  n = 0;
  prev = suffix;
}


{
  s [n ++] = $0;
}


END {
  if (prev == "-001")
    Pr();
}


function Pr ()
{
  for (i = 0; i < n; i ++)
    print s [i];
}


