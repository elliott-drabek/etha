#!/bin/awk -f
# Usage:  exon-ends.awk
#   Read a sorted file of putative exon end positions and output
#   the best one from each cluster.  Input format is, on each line:
#     sequence id
#     putative end base position (counting from 1)
#     direction of exon, either "f" or "r"
#     sequence length
#     end precision; this is the number of bases imputed to get the
#       end position.  A lower number here indicates a more accurate
#       end position.
#   Lines must be sorted -k1,1b -k3,3b -k2nb


BEGIN {
  Sep = 500;
    # consecutive putative ends at least this close are in the same cluster
}


{
  if ($2 < 1 || $4 < $2)
    next;  # skip entries off the end of the sequence

  if ($1 == id && $3 == dir && $2 - end <= Sep)
    {  # same cluster keep if better than what's saved
      if ($5 < prec)
        {
          best_line = $0;
          best_end = $2;
          prec = $5;
        }
      end = $2;
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
  end = best_end = $2;
  dir = $3;
  prec = $5;
}


