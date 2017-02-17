#!/bin/awk -f
# Usage:  exon-starts.awk
#   Read a sorted file of putative exon start positions and output
#   the best one from each cluster.  Input format is, on each line:
#     sequence id
#     putative start position
#     direction from start, either "f" or "r"
#     sequence length
#     start precision; this is the number of bases imputed to get the
#       start position.  A lower number here indicates a more accurate
#       start position.
#   Lines must be sorted -k1,1b -k3,3b -k2nb


BEGIN {
  Sep = 500;
    # consecutive putative starts at least this close are in the same cluster
}


{
  if ($2 < 1 || $4 < $2)
    next;  # skip entries off the end of the sequence

  if ($1 == id && $3 == dir && $2 - start <= Sep)
    {  # same cluster keep if better than what's saved
      if ($5 < prec)
        {
          best_line = $0;
          best_start = $2;
          prec = $5;
        }
      start = $2;
    }
  else
    {
      Pr();
    }
}


END {
  Pr();
}


function Pr ()
{
  if (best_line != "")
    printf "%s\n", best_line;
  best_line = $0;
  id = $1;
  start = best_start = $2;
  dir = $3;
  prec = $5;
}


