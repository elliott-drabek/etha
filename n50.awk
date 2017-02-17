#!/bin/awk -f
# Usage:  n50.awk [-v target=<n>]
# Reads a file from stdin and computes the N50 value
# Input MUST be in descending sorted order
# If the value of variable  target  is set (by the -v
# option) then that value will be used to compute N50
# instead of the sum of the inputs

BEGIN   {
  n = 0;
}


{
  a [n] = $1;
  n ++;
  s += $1;
}


END     {
  if  (target == 0)
    target = s / 2.0;
  else
    target /= 2.0;

  printed = 0;
  for  (i = 0;  i < n && t < 2.0 * target;  i ++)
    {
      t += a [i];
      if  (target <= t && ! printed)
        {
          printf "total n = %d  i = %d  s = %.0f  target = %.0f  n50 = %.0f\n",
            n, i + 1, s, 2.0 * target, a [i];
          printf "  sum = %.0f  prior = %.0f\n",
            t, a [i - 1];
          printed = 1;
        }
    }

  j = int (i / 2);
  if  (i % 2 == 1)
    median = a [j];
  else
    median = (a [j] + a [j - 1]) / 2.0;
  printf "used n = %d  mean = %.0f  median = %.0f  max = %d  min = %d\n",
    i, t / i, median, a [0], a [i - 1];
} 

