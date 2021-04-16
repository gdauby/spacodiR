



#' @title Masuring spatial and phylogenetic structuring of diversity in communities
#'
#' @description Considering species-, phylogenetic-, or trait-diversities, this function measures diversity structuring of community samples.
#'
#' @param sp.plot a community dataset in spacodiR format (see as.spacodi)
#' @param phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)
#' @param sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
#' @param all.together Boolean; whether to treat all traits together or separately
#' @param prune Boolean whether to dynamically prune datasets if mismatches occur
#' @param pairwise Boolean whether to return pairwise diversity measures amongst all plots
#' @param ... Additional arguments to be passed to match.spacodi.data
#' 
#' @author Olivier Hardy, Timothy Paine, and Jonathan Eastman
#' 
#' @details 
#' spacodi.calc requires a community dataset (species-by-plots matrix; sp.plot) of numerical abundance, 
#' relative abundance, or presence | absence for plots. spacodi.calc returns statistics of diversity partitioning of plots, 
#' considering species diversity and, if additional information is provided, either trait or phylogenetic diversities among plots. 
#' If phy=NULL and sp.traits=NULL, a measure of partitioning for species diversity will be returned.
#' In treating each pair of plots as a community unto its own, pairwise=TRUE will return estimates for diversity structuring 
#' for all pairwise combinations of plots.
#' 
#' If a phylogeny or trait dataset is supplied with species that are not present in the community dataset (i.e., sp.plot) 
#' or vice versa, the user has the option to dynamically prune these datasets to match (prune=TRUE). 
#' 
#' If prune=FALSE and dataset mismatches occur, the function will inevitably return NaN where plots have fewer than two distinct species sampled.
#' 
#' @return 
#' A named list of at least one element (Ist) is returned. The size of the returned list is wholly dependent upon given arguments.
#' 
#' SPECIES DIVERSITY STRUCTURING
#' 
#' \itemize{
#'   \item Ist: a measure of local species identity excess between individuals, expressing species turnover. 
#'   It is a form of spatial partition of Gini-Simpson diversity (equivalent to Fst in population genetics). 
#'   Ist considers only abundances (or presences) in the species-by-plots matrix.
#' }
#' 
#' PHYLOGENETIC DIVERSITY STRUCTURING
#' 
#' \itemize{
#'   \item Pst: a measure of local phyletic proximity excess between individuals, expressing species + phylogenetic turnover. 
#'   It is a form of spatial partition of Rao's quadratic entropy (equivalent to Nst in population genetics). 
#'   Tst is the analogue for trait data, estimating the spatial partitioning of mean trait-divergence between individuals.
#' }
#' 
#' @references 
#' HARDY OJ and B SENTERRE. 2007. Characterizing the phylogenetic structure of communities by an additive partitioning of phylogenetic diversity. Journal of Ecology 95:493-506.
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' HARDY OJ and L JOST. 2008. Interpreting and estimating measures of community phylogenetic structuring. Journal of Ecology 96:849-852.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' 
#' # community diversity statistics of Hardy and Senterre (2007): tree-based
#' spacodi.calc(sp.plot = spl, phy = phy)
#' 
#' # community diversity statistics: trait-based with pairwise comparisons
#' spacodi.calc(sp.plot = spl, phy = phy, pairwise=TRUE)
#' 
#' # community diversity for a pair of traits
#' spacodi.calc(sp.plot = spl, sp.traits = trt, all.together=TRUE)
#' 
#' 
#' # community diversity for a pair of traits, each singly
#' spacodi.calc(sp.plot = spl, sp.traits = trt, all.together=FALSE)
#' 
#' # Ist: using abundance data only
#' spacodi.calc(sp.plot = spl)
#' 
#' # calculations with missing taxa between tree and sp.plot
#' # excluding the last five species in sp.plot, 
#' spacodi.calc(sp.plot = spl[1:15,], phy = phy, prune=TRUE)
#' 
#' # as before but with 'manual' pruning of the datasets
#' match.spacodi.data(sp.plot=spl[1:15,],phy=phy) -> prn.data
#' spacodi.calc(sp.plot=prn.data$sp.plot, phy=prn.data$sp.tree)
#' prn.data$sp.plot
#' prn.data$sp.tree
#' 
#' @importFrom stats dist
#' 
#' @export
spacodi.calc <-
  function(sp.plot,
           phy = NULL,
           sp.traits = NULL,
           all.together = TRUE,
           prune = TRUE,
           pairwise = FALSE,
           ...) {
    
    if (!missing(phy) &&
        !missing(sp.traits))
      stop("Please supply either a phylogeny or a trait matrix, but not both.")
    
    
    # determine type of abundance
    stripped = unname(unlist(c(sp.plot)))
    if (all(!is.na(match(stripped, c(0, 1))))) {
      abundt = 0		# presence|absence
    } else if (all(!is.na(match(range(stripped), c(0, 1)))) &&
               length(unique(stripped)) > 2) {
      abundt = 1		# relative abundance
    } else {
      abundt = 2		# abundance (n.individuals)
    }
    
    # iter is 1 unless more than a single trait is being used for separate analyses
    iter <- 1
    traits.tmp <- distmat <- NULL
    
    # INTERPRET whether to compute trait diversity or phylogenetic diversity & prepare data
    if (!missing(phy)) {
      # distmat is phylogenetic: Pst
      sp.data <- match.spacodi.data(sp.plot = sp.plot,
                                   phy = phy,
                                   prune = prune,
                                   ...)
      sp.plot <- sp.data$sp.plot
      phy <- sp.data$sp.tree
      if (check.distmat(phy))
        distmat <- phy
      else
        distmat <- cophen(phy)
      
    } else if (missing(sp.traits) & missing(phy)) {
      # distmat is null: Ist
      distmat <- matrix(1, ncol = nrow(sp.plot), nrow = nrow(sp.plot))
      
    } else if (!missing(sp.traits)) {
      # distmat is trait-based: Tst
      if (class(sp.traits) == "dist")
        sp.traits <- as.matrix(sp.traits)
      if (ncol(sp.traits) == 1)
        all.together <- TRUE
      if (all(is.null(names(sp.traits))))
        names(sp.traits) <- paste("trt", seq(1:ncol(sp.traits)), sep = "")
      if (all.together == TRUE | check.distmat(sp.traits)) {
        sp.data <- match.spacodi.data(sp.plot = sp.plot,
                                     sp.traits = sp.traits,
                                     prune = prune,
                                     ...)
        sp.plot <- sp.data$sp.plot
        sp.traits <- sp.data$sp.traits
        if (check.distmat(sp.traits))
          distmat <- sp.traits
        else
          distmat <- as.matrix(stats::dist(sp.traits))
      } else {
        iter <- ncol(sp.traits)
        traits.tmp <- lapply(1:ncol(sp.traits), function(x) {
          trait.tt <- data.frame(sp.traits[, x])
          row.names(trait.tt) <- row.names(sp.traits)
          sp.data <-
            match.spacodi.data(sp.plot = sp.plot,
                               sp.traits = trait.tt,
                               prune = prune,
                               ...)
          distmat <- as.matrix(stats::dist(sp.data$sp.traits))
          return(list(distmat = distmat, sp.plot = sp.data$sp.plot))
        })
      }
    } else if (is.null(distmat)) {
      stop("Cannot decipher input object(s).")
    }
    
    if (is.null(names(sp.plot)))
      pnames <- paste("plt", seq(1, ncol(sp.plot)), sep = ".")
    else
      pnames <- names(sp.plot)
    
    
    # PREPARE output
    gen.out <- list()
    prw.out <- list()
    
    for (tt in 1:iter) {
      if (is.null(traits.tmp)) {
        sp.plot <-	as.matrix(sp.plot)
      } else {
        sp.plot <-	as.matrix(traits.tmp[[tt]]$sp.plot)
        distmat <-	as.matrix(traits.tmp[[tt]]$distmat)
      }
      
      diag(distmat) <- 0
      np <- ncol(sp.plot)
      out <- NA
      out <- .C(
        "spacodi",
        np = as.integer(np),
        ns = as.integer(nrow(sp.plot)),
        sp.plot = as.double(as.vector(sp.plot)),
        distmat = as.double(as.vector(as.matrix(distmat))),
        abundtype = as.integer(abundt),
        Ndclass = as.integer(0),
        dclass = as.double(c(0, 0)),
        Ist = 0,
        Pst = 0,
        Bst = 0,
        PIst = 0,
        pairwiseIst = as.double(as.vector(matrix(0, np, np))),
        pairwisePst = as.double(as.vector(matrix(0, np, np))),
        pairwiseBst = as.double(as.vector(matrix(0, np, np))),
        pairwisePIst = as.double(as.vector(matrix(0, np, np))),
        PACKAGE = "spacodiR"
      )
      
      # compile results for phylogenetic turnover
      if (missing(phy) & missing(sp.traits)) {
        r.out <- as.numeric(out[8])
        names(r.out) <- "Ist"
        gen.out[[tt]] <- as.data.frame(t(r.out))
        if (pairwise) {
          prw.out[[tt]] <- list(rematrix(unlist(out[12]), pnames))
          names(prw.out[[tt]]) <- paste("pairwise", "Ist", sep = ".")
        }
      } else if (!missing(phy)) {
        r.out <- as.numeric(c(out[8:11]))
        names(r.out) <- c("Ist", "Pst", "Bst", "PIst")
        gen.out[[tt]] <- as.data.frame(t(r.out))
        if (pairwise) {
          prw.out[[tt]] <- lapply(out[12:15], function(x)
            rematrix(x, pnames))
          names(prw.out[[tt]]) <-
            paste("pairwise", c("Ist", "Pst", "Bst", "PIst"), sep = ".")
        }
      } else if (!missing(sp.traits)) {
        r.out <- as.numeric(c(out[8:11]))
        names(r.out) <- c("Ist", "Tst", "Ust", "TAUst")
        gen.out[[tt]] <- as.data.frame(t(r.out))
        if (pairwise) {
          prw.out[[tt]] <- lapply(out[12:15], function(x)
            rematrix(x, pnames))
          names(prw.out[[tt]]) <-
            paste("pairwise", c("Ist", "Tst", "Ust", "TAUst"), sep = ".")
        }
      }
    }
    if (iter > 1 & !all.together) {
      RES = lapply(1:iter, function(x) {
        if (pairwise)
          return(c(gen.out[[x]], prw.out[[x]]))
        else
          return(c(gen.out[[x]]))
      })
      names(RES) <- names(sp.traits)
    } else {
      prw.out <- unlist(prw.out, recursive = FALSE)
      gen.out <- unlist(gen.out, recursive = FALSE)
      if (pairwise)
        RES <- c(gen.out, prw.out)
      else
        RES <- gen.out
    }
    
    return(unlist(
      list(
        RES,
        sp.plot = list(sp.plot),
        sp.tree = list(phy),
        sp.traits = list(sp.traits)
      ),
      recursive = FALSE
    ))
    
  }

		
#' @title Internal function
#'
#' @description Internal function
#' 
#' @param arr vector
#' @param names character vector of species names
#' 
#' @author JM Eastman
#' 
#' @return a phylo object
rematrix <- function(arr, names) {
  n = sqrt(length(arr))
  counter = 1
  M = matrix(0, n, n)
  for (i in 1:n) {
    j = i + 1
    while (j <= n) {
      M[i, j] = arr[counter]
      counter = counter + 1
      j = j + 1
    }
  }
  M[lower.tri(M)] = t(M)[lower.tri(t(M))]
  dimnames(M)[[1]] <- dimnames(M)[[2]] <- names
  return(M)
}