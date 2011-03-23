This directory includes a small example of how METAL can be used to meta-analyze study results. 

The metal script file is called metal.txt. To run metal with this input file, issue one of the
following commands:

  # Using input redirection
  metal < metal.txt

  # Using a command line script parameter
  metal metal.txt

  # Using metal in interactive mode
  metal [ENTER]
  SOURCE metal.txt
  QUIT

The example is drawn from a recent study of the genetics of glucose levels (Chen et al, Journal of
Clinical Investigation, 2008; Prokopenko et al, Nature Genetics, 2009). Input files summarizing
results for each of three studies (in a subset of interesting regions) are called:

  DGI_three_regions.txt 
  MAGIC_FUSION_Results.txt.gz 
  magic_SARDINIA.tbl

Each of the files uses a slightly different format to report results and this is accomodated in 
the metal.txt script.


