#!/bin/awk -f
# Usage:  extract-fasta-bytag.awk  <list-file>
#   Read a multifasta file and output records in it that
#   have their tag in <list-file>.  The tag is the rest
#   of field one after the '>'

BEGIN   {
         if  (ARGC < 2)
             Usage_Exit();

         fp = ARGV [1];
         delete ARGV [1];
         while  ((getline < fp) > 0)
           tag [$1] = 1;
        }

        {
         if  (substr ($1, 1, 1) == ">")
             s = tag [substr ($1, 2)];

         if  (s == 1)
             print;
        }


function  Usage_Exit  ()
  {
   print "# Usage:  extract-fasta-bytag.awk  <list-file>";
   print "#   Read a multifasta file and output records in it that";
   print "#   have their tag in <list-file>.  The tag is the rest";
   print "#   of field one after the '>'";

   exit;
  }


