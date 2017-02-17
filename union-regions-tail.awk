#!/bin/awk -f
# Usage:  union-regions-tail.awk
#   Read sorted list of regions from stdin and output the tail of
#   joined overlapping regions, where the tail is at least Tail_Len
#   if possible and starts at one of the input region starts.
#   Input format is:
#     id   lo   hi   dir
#   where dir is 'f' or 'r' to indicate forward or reverse strand.
#   Input is assumed sorted -k1,1 -k2n -k3nr.
#   Output is the same format as the input.


BEGIN {
  Tail_Len = 500;
  p = "";
}


{
  if ($1 != p || $4 != dir || hi < $2)
    Pr();
  else
    {
      if (hi < $3)
        hi = $3;
      if (Tail_Len <= 1 + hi - $2)
        lo = $2;
    }
}


END {
  Pr();
}


function Pr ()
{
  if (p != "")
    printf "%s %7d %7d  %s\n", p, lo, hi, dir;
  p = $1;
  lo = $2;
  hi = $3;
  dir = $4;
}
