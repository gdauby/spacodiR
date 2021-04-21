

#' @title Adding species to congeneric ones in a phylogenetic tree
#'
#' @description Adding species on a phylogenetic tree when the genus is present
#'
#' @param specieslist a community dataset in spacodiR format (see as.spacodi)
#' @param phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)
#' @param threshold_upper_branch a numeric value indicating what is the maximal branch length for graftin the species, default is 20
#' 
#' @author Olivier Hardy, Gilles Dauby
#' 
#' @details 
#' The initial tree must have tip labels in the form of "genus_species")
#' A list of species to add must be given (specieslist)
#' Each species is added in turn as sister to a randomly chosen congeneric species
#' adding a new node with branch length to tips = 0.9 the initial branch length
#' but with a maximal branch length of 20. => to improve, e.g. considering the distribution
#' of distance between congeneric spp in the original tree
#' 
#' 
#' 
#' @return 
#' A named list
#' \itemize{
#'   \item phylo:  the completed tree 
#'   \item species_already_present: species already present in the initial tree
#'   \item non_added_species species not added because of lack of congeneric species in original tree
#' }
#'  
#' 
#' @examples 
#' 
#' @import ape
#' 
#' @export
add_species_to_genera <- function(specieslist, phy, threshold_upper_branch = 20)
{
  #extract genus for each species in phylogeny
  initial_species <- phy$tip.label
  generatree <-
    sapply(strsplit(phy$tip.label, split = "_"), '[', 1) #extract genus name from species in tree
  generalist <-
    sapply(strsplit(specieslist, split = "_"), '[', 1) #extract genus name from species list
  
  species_already_present <- specieslist %in% phy$tip.label #False/True vector
  genera_already_present <- generalist %in% generatree #False/True vector
  species_to_add <- specieslist[genera_already_present & !(species_already_present)]  #list of species to add (species non existing in the tree but from genus already in the tree)
  genus_to_add <- generalist[genera_already_present & !(species_already_present)] #idem but with genus name of each species
  non_added_species <- specieslist[!genera_already_present] #list of species that cannot be added
  species_already_present <- specieslist[species_already_present] #list of species already
  
  #write done numbers
  cat(paste(length(species_already_present),"species already present,", length(non_added_species),"species without matching genus (not added),", length(species_to_add),"species to add now to the tree: "))
  
  
  for (splist in 1:length(species_to_add)) { #loop over species in the list
    #    order=c(1:length(tree$tip.label))
    newsp <- species_to_add[splist]
    genus <- genus_to_add[splist]
    
    tips <-  which(generatree == genus) #list the indices of species in the tree of given genus
    if (length(tips) > 1) {
      tip <-
        sample(tips, 1)      
    } else {
      tip <- tips #chose one random congeneric species from the phylogeny
    }
    
    postipinedge <-
      which(phy$edge[, 2] == tip)[1] #identify position of chosen tip in matrix edge (if there are several take only the first one => to improve by warning that not all species names are different)
    
    upnodeoftip <-
      phy$edge[postipinedge, 1]  #identify upper node of chosen tip
    lentip <-
      phy$edge.length[postipinedge] #length of branch to chosen tip
    newlentip <-
      lentip * 0.9 #add the new species at a dist of 0.9 of the upper noderunif(1) #chose a uniform branch length to branch the new species
    if (newlentip > threshold_upper_branch)
      newlentip <- threshold_upper_branch #put a max diverg time between congeners (to improve by computing the mean distance between congeners in original tree)
    
    tr2 <-
      phy  #create new phylotree with a new node for the added species
    tr2$tip.label <-
      c(phy$tip.label, newsp)  #copy tips and add the new one at end
    Ntips <-
      1 + length(phy$tip.label) #Nb of tips in new tree (and number of new tip to add)
    generatree <-
      c(generatree, genus) #update the list of genera in the tree
    tempedge <- phy$edge #create a temporary table for edges
    tempedge[, 1] <-
      tempedge[, 1] + 1 #add one to each node number as a new species has been added (1st column)
    tempedge[tempedge[, 2] > Ntips, 2] <-
      tempedge[tempedge[, 2] > Ntips, 2] + 1 #same for 2nd column but without changing tip numbers
    tr2$Nnode <- as.integer(phy$Nnode + 1) #add a node
    newnode <-
      Ntips + tr2$Nnode #number of new node (should be the last number of the parent tree + 2)
    
    edges_to_add <-
      rbind(
            c(upnodeoftip + 1, newnode),
            c(newnode, tip),
            c(newnode, Ntips)) #replace the edge of the tip by 3 edges (the one going from the new node to the ancestral node plus the tip and new tip to the new node)
    
    if(postipinedge > 1) {
      
      edges_to_add <- 
        rbind(tempedge[1:(postipinedge - 1), ],
              edges_to_add)
      
    }
    
    tr2$edge <-
      edges_to_add
    
    edges_length_to_add <- 
      c(lentip - newlentip, newlentip, newlentip)
    
    if(postipinedge > 1) {
      
      edges_length_to_add <- 
        c(phy$edge.length[1:(postipinedge - 1)], edges_length_to_add)
      
    }
      
    # tr2$edge.length <-
    #   c(phy$edge.length[1:(postipinedge - 1)], lentip - newlentip, newlentip, newlentip) #copy tip branches lengths
    
    tr2$edge.length <- edges_length_to_add
    
    if (postipinedge < length(tempedge[, 1])) {
      
      tr2$edge <-
        rbind(tr2$edge, tempedge[(postipinedge + 1):length(tempedge[, 1]), ]) #add rest of the edge
      tr2$edge.length <-
        c(tr2$edge.length, phy$edge.length[(postipinedge + 1):length(phy$edge[, 1])]) #copy tip branches lengths
      
    }
    
    phy <- tr2
    #   if(splist<1000) if(splist/100==round(splist/100)) cat(splist,"")
    if (splist / 100 == round(splist / 100))
      cat(splist, "")
  }
  output <- list(phy, non_added_species, species_already_present)
  names(output) <- c("phylo","non_added_species","species_already_present")
  return(output)
}
