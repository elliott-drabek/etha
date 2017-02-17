#!/bin/awk -f
# Usage:  make-unitig-seq.awk  <fa>
#   Read sequences from fasta file <fa>.  Then read unitig information
#   from stdin (in the format produced by program unitig) that references
#   those sequences and output the corresponding unitig sequences.


BEGIN {
  if  (ARGC < 2)
    Usage_Exit();

  fp = ARGV [1];
  delete ARGV [1];

  s = "";
  while  ((getline < fp) > 0)
    {
      if (match ($0, /^>/))
        {
          if (s != "")
            seq [tag] = s;
          s = "";
          tag = substr ($1,2);
        }
      else
        s = s $1;
    }
  if (s != "")
    seq [tag] = s;

  s = "";
  in_unitig = 0;
}


/^Unitig/ {
  if (in_unitig)
    Output();

  id = substr ($2, 2);
  in_unitig = 1;
  next;
}


{
  if (in_unitig)
    {
      if (NF == 0)
        {
          Output();
          in_unitig = 0;
        }
      else if (NF == 5)
        {
          if (seq [$5] == "")
            {
              printf "ERROR:  sequence %s not found\n", $5;
              exit -1;
            }
          else if (length (seq [$5]) != $3)
            {
              printf "ERROR:  end position %d != sequence length %d\n",
                $3, length (seq [$5]);
              printf "  line:  %s\n", $0;
              exit -1;
            }
          s = s substr (seq [$5], $2);
        }
    }
  
}


END {
  if (in_unitig)
    Output();
}


function Output ()
  # Output the current unitig
{
  printf ">uni%04d\n", id;

  n = length (s);
  for (i = 1; i <= n; i += 60)
    print substr (s, i, 60);

  s = "";
  return;
}


function Usage_Exit ()
{
  print "# Usage:  make-unitig-seq.awk  <fa>";
  print "#   Read sequences from fasta file <fa>.  Then read unitig information";
  print "#   from stdin (in the format produced by program unitig) that references";
  print "#   those sequences and output the corresponding unitig sequences.";

  exit;
}
