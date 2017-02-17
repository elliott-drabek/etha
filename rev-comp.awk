#!/bin/awk -f

# Reverse complement the fasta input string
        {
         if  (substr($1,1,1) == ">")
             {
              Rev_Print(tag, s);
              tag = $0 "  comp";
              s = "";
             }
           else
             {
              gsub (/[ \t]/, "");
              s = s $0;
             }
        }
END     {
         Rev_Print(tag, s);
        }


function  Rev_Print  (tag, s)
  {
   len = length (s);

   if  (len == 0)
       return;
   print tag;
   rev = "";
   for  (i = len;  i > 0;  i --)
     rev = rev Complement(substr (s, i, 1));

   for  (i = 1;  i <= len;  i += 60)
     {
      if  (1 + len - i > 60)
          j = 60;
        else
          j = 1 + len - i;
      print substr (rev, i, j);
     }

   return;
  }


function  Complement  (s)
  {
   if  (s == "A")
       return  "T";
   if  (s == "C")
       return  "G";
   if  (s == "G")
       return  "C";
   if  (s == "T")
       return  "A";

   if  (s == "R")    # A or G
       return  "Y";
   if  (s == "Y")    # C or T
       return  "R";
   if  (s == "S")    # G or C
       return  "S";
   if  (s == "W")    # A or T
       return  "W";
   if  (s == "K")    # G or T
       return  "M";
   if  (s == "M")    # A or C
       return  "K";

   if  (s == "B")    # C or G or T
       return  "V";
   if  (s == "D")    # A or G or T
       return  "H";
   if  (s == "H")    # A or C or T
       return  "D";
   if  (s == "V")    # A or C or G
       return  "B";

   if  (s == "N")
       return  "N";

   if  (s == "a")
       return  "t";
   if  (s == "c")
       return  "g";
   if  (s == "g")
       return  "c";
   if  (s == "t")
       return  "a";
   if  (s == " ")
       return  "";
   if  (s == "r")    # a or g
       return  "y";
   if  (s == "y")    # c or t
       return  "r";
   if  (s == "s")    # g or c
       return  "s";
   if  (s == "w")    # a or t
       return  "w";
   if  (s == "k")    # g or t
       return  "m";
   if  (s == "m")    # a or c
       return  "k";

   if  (s == "b")    # c or g or t
       return  "v";
   if  (s == "d")    # a or g or t
       return  "h";
   if  (s == "h")    # a or c or t
       return  "d";
   if  (s == "v")    # a or c or g
       return  "b";

   return  "n";
  }



