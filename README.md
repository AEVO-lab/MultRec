------------------------------------------------------------------
MULTREC - Multi-reconciliation program 
------------------------------------------------------------------
Multrec takes as input a species tree S, a set of gene trees, a duplication cost, a loss cost and a parameter duplication height h.  The output is a mapping of the gene tree nodes to S that minimizes the segmental reconciliation cost, assuming that there exists such a mapping that has duplications sum-of-heights at most h.  If loss cost >= dup cost, the LCA mapping is returned.
The leaves of the gene trees must map to the leaves of S.  The gene tree leaves are assumed to have the format [species_name]__[gene_name], for example HUMAN_BRCA2 indicates that the gene is mapped to the leaf of S names HUMAN.  The gene/species separator can be changed with the -spsep argument, and the position of the species name in the gene name with the -spindex argument, indexed at 0.  
If your genes are name e.g. GENENAME_SPECIESNAME_OTHERSTUFF, you can set -spsep \_\ -spindex 1

The format of the output is 5 lines, as follows
COST=[total cost of mapping]
DUPHEIGHT=[sum of duplication heights]
NBLOSSES=[number of losses]
SPECIESTREE=[species tree newick, with internal nodes labeled
GENETREES=[all gene tree newick separated by ;.  Internal nodes are labeled by mapping]

If no solution is found (when h is too small), then the output is simply
NO SOLUTION FOUND

Here is an example execution
./Multrec -g "((A__1, C__1),B__1);((A__2, B__2),B__3);" -s "((A,B),(C,D));" -d 3

Required arguments:
At least one of -g or -gf must be specified, and at least one of -s or -sf must be specified.
-g   [g1;g2;...;gk]   Here g1,g2,...,gk are gene trees
                      represented in Newick format.  
                      The gene trees are separated by the ; symbol.	
-gf  [file]           file is the name of a file containing the list 
                      of gene trees, all in Newick format and separated 
                      by a ; symbol in the file.
-s   [newick]         The species tree in Newick format.
-sf  [file]           Name of the file containing species tree Newick.

Optional arguments:
--help                Print this help message.
-d   [double]         The cost for one height of duplication.  Default=3
-l   [double]         The cost for one loss.  Default=1
-h   [int]            Maximum allowed duplication sum-of-heights.  Default=20
-o   [file]           Output file.  Default=output to console
-spsep   [string]     Gene/species separator in the gene names.  Default=__
-spindex [int]        Position of the species in the gene names, after 
                      being split by the gene/species separator.  Default=0
--test                Launches a series of unit tests.  This includes small fixed 
                      examples with known outputs to expect, and larger random trees 
                      to see if the program terminates in an OK status on more complicated
                      datasets.
