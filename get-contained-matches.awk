#!/bin/awk -f
# Usage:  get-contained-matches.awk
#   Read output from show-coords -cHlT and output references that match
#   queries end-to-end or contain queries as substrings.


{
  if ($10 == 100.0 && $11 == 100.0 && 99.0 <= $7)
    printf "%-15s %-8s  %-8s\n", $12, "equals", $13;
  else if ($3 == 1 && $4 == $9 && 99.0 <= $7)
    printf "%-15s %-8s  %-8s\n", $12, "contains", $13;
}


