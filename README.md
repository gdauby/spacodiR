

This repository contains the R package spacodiR which was removed from CRAN.

Some functions were removed as they were not used or necessary for the key function `spacodi.calc` used to compute various metrics of phylogenetic/functional community structures.


All functions were written by JE Eastman, Tim CE Paine and OJ Hardy.

Computations of various statistics rely on C code also found in an executable software that can be found in the following link:

https://ebe.ulb.ac.be/ebe/SPACoDi.html  

The package can be installed by the following code

`install.packages("devtools")`
`devtools::install_github("gdauby/spacodiR")`


See literature below for further details:

A formal description of the methods  
HARDY OJ and B SENTERRE. 2007. Characterizing the phylogenetic structure of communities by an additive partitioning of phylogenetic diversity. Journal of Ecology 95:493-506. 

Statistical properties of the approach and permuation tests  
HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.

Additional details on the interpretation  
HARDY OJ and L JOST. 2008. Interpreting and estimating measures of community phylogenetic structuring. Journal of Ecology 96:849-852.

A presentation of the original package  
EASTMAN JE, Paine TCE, Hardy OJ 2011 spacodiR: structuring of phylogenetic diversity in ecological communities. Bioinformatics, Volume 27, Issue 17, 1 September 2011, Pages 2437–2438,


