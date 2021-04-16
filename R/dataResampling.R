#' @title Rndomizing a community phylogenetics matrix: '1a' of Hardy (2008)
#'
#' @description  Used for resampling data within a community dataset
#'
#' @param obj a community dataset in spacodiR format (see as.spacodi)
#' @param abund.class.ratio numeric vector of length 1 defining abundance classes
#' 
#' @author Timothy Paine and Jonathan Eastman
#' 
#' @details 
#' A resampling procedure for a species-by-plots matrix, where species are shuffled within abundance classes. 
#' Species are grouped into distinct abundance classes characterized by a fixed ratio: abund.class.ratio = maximal abundance / minimal abundance. 
#' For instance, if abund.class.ratio = 4, the limits between abundance classes could be 1, 4, 16, ... . Species are randomly permuted 
#' within each class, which maintains most of the abundance phylogenetic structure originally present in a dataset.
#' 
#' @return A shuffled dataset
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' # shuffle dataset
#' resamp.1a(obj=spl, abund.class.ratio=3) 
#' 
#' @importFrom stats runif
#' 
#' 
#' @export
resamp.1a <-
  function(obj, abund.class.ratio = 4) {
    if (abund.class.ratio <= 1)
      stop("Supplied abundance-class ratio did not appear sensible: choose a value greater than 1.")
    orig <- obj
    abund <- rowSums(obj)
    n.spp   <- length(abund)
    
    aa <- abund.class.ratio
    
    while (1) {
      classes <- aa ^ (0:ifelse(aa < 2, n.spp * (1 / (aa - 1)), n.spp))
      if (length(classes) >= n.spp) {
        classes <- unique(round(runif(1) * classes))
        break()
      }
    }
    
    class <- rep(NA, n.spp)
    for (i in 1:length(classes)) {
      class[abund > classes[i]] <- i
    }
    if (any(is.na(class)))
      class[which(is.na(class))] = 1
    new_name <- rep(NA, n.spp)
    for (i in unique(class)) {
      new_name[class == i] <- sample(rownames(obj[class == i, ]))
    }
    row.names(obj) <- new_name
    
    return(obj[order(match(row.names(obj), row.names(orig))), ])
  }


#' @title randomizing a community phylogenetics matrix: '1s' of Hardy (2008)
#'
#' @description Used for resampling data within a community dataset
#'
#' @param obj a community dataset in spacodiR format (see as.spacodi)
#' 
#' @author Timothy Paine and Jonathan Eastman
#' 
#' @details 
#' A resampling procedure for a species-by-plots matrix, where observed abundances are shuffled across species and plots
#' 
#' @return A shuffled dataset
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' 
#' # shuffle dataset
#' resamp.1s(obj=spl) 
#' 
#' @export
resamp.1s <-
function(obj) {
	orig=obj
	row.names(obj) <- sample(row.names(obj))
	return(obj[order(match(row.names(obj),row.names(orig))),])
}


#' @title randomizing a community phylogenetics matrix: '2s' of Hardy (2008)
#'
#' @description Used for resampling data within a community dataset
#'
#' @param obj a community dataset in spacodiR format (see as.spacodi)
#' 
#' @author Timothy Paine and Jonathan Eastman
#' 
#' @details 
#' A resampling procedure for a species-by-plots matrix, where observed abundances are shuffled across species but within plots.
#' 
#' @return A shuffled dataset
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' 
#' # shuffle dataset
#' resamp.2s(obj=spl) 
#' 
#' @export
resamp.2s <-
function(obj) {
	for(nn in 1:ncol(obj)){
		obj[,nn]=sample(obj[,nn])
	}
	return(obj)
}


#' @title randomizing a community phylogenetics matrix: '2x' of Hardy (2008)
#'
#' @description Used for resampling data within a community dataset
#'
#' @param obj a community dataset in spacodiR format (see as.spacodi)
#' @param level numeric vector of length one giving a proportion specifying the extent of data shuffling
#' 
#' @author Timothy Paine and Jonathan Eastman
#' 
#' @details 
#' A resampling procedure for a species-by-plots matrix, based on Gotelli swapping. 
#' Shuffles abundances within a pair of plots and for a pair of species. 
#' The level defines the degree of sampling, with larger values dictating a higher level of reshuffling. 
#' For instance, if level = 0.4 and the dataset involves 5 species and 10 plots, a total of 20 (0.4x5x10) Gotelli swaps are performed
#' 
#' @return A shuffled dataset
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' GOTELLI NJ. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606-2621.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' 
#' # shuffle dataset
#' resamp.2x(obj=spl, level=0.2) 
#' 
#' @export
resamp.2x <-
  function(obj, level = 0.1) {
    swaps = round(level * ncol(obj) * nrow(obj))
    orig = obj
    for (swap in 1:swaps) {
      rcol = sample(1:ncol(obj), 2)
      rspp = sample(1:nrow(obj), 2)
      obj[rspp[1], rcol[1]] = orig[rspp[2], rcol[1]]
      obj[rspp[2], rcol[1]] = orig[rspp[1], rcol[1]]
      obj[rspp[1], rcol[2]] = orig[rspp[2], rcol[2]]
      obj[rspp[2], rcol[2]] = orig[rspp[1], rcol[2]]
      orig = obj
    }
    if (sum(obj) != sum(orig))
      warning("A poor result is likely.")
    return(obj)
  }


#' @title randomizing a community phylogenetics matrix: '3i' of Hardy (2008)
#'
#' @description Used for resampling data within a community dataset
#'
#' @param obj a community dataset in spacodiR format (see as.spacodi)
#' 
#' @author Jonathan Eastman
#' 
#' @details 
#' A resampling procedure for a species-by-plots matrix, where observed abundances are shuffled within species and among plots.
#' 
#' @return A shuffled dataset
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' 
#' # shuffle dataset
#' resamp.3i(obj=spl)
#' 
#' @export
resamp.3i <-
  function(obj) {
    for (ss in 1:nrow(obj)) {
      obj[ss, ] = sample(obj[ss, ])
    }
    return(obj)
  }


#' @title randomizing a community phylogenetics matrix: '3t' of Hardy (2008)
#'
#' @description Used for resampling data within a community dataset
#'
#' @param obj a community dataset in spacodiR format (see as.spacodi)
#' @param dmat an optional dataframe of distances between plots; row names and column names should be identical
#' 
#' @author Jonathan Eastman
#' 
#' @details 
#' A resampling procedure for a species-by-plots matrix, where observed abundances within species are shuffled to adjacent plots. 
#' This procedure thus assumes meaningful arrangement of plots in space. If a distance matrix is supplied, 
#' the likelihood of shuffling to a particular plot is proportional to the distance between the plots.
#' 
#' @return A shuffled dataset
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' 
#' # define a distance matrix
#' foo <- matrix(runif((ncol(spl)->ss)^2,0,100),ss,ss)
#' foo[upper.tri(foo)] <- foo[lower.tri(foo)]
#' diag(foo) <- 0
#' dmat <- as.data.frame(foo)
#' row.names(dmat) <- names(spl)
#' names(dmat) <- row.names(dmat)
#' 
#' 
#' # shuffle dataset
#' resamp.3t(obj=spl, dmat=dmat) 
#' spl ## comparison with original
#' 
#' @importFrom stats runif
#' 
#' @export
resamp.3t <-
function(obj, dmat=NULL) {
	if(is.null(dmat))flag=TRUE else flag=FALSE
	names.orig=names(obj)
	names(obj)=seq(1:ncol(obj))
	if(is.null(dmat)) {
		dmat=as.data.frame(matrix(0,ncol(obj),ncol(obj)))
		names(dmat)=names.orig
		row.names(dmat)=names.orig
	}
	dmat=as.data.frame(as.matrix(dmat))
	if(!all(names(dmat)%in%names.orig) || ncol(obj)!=ncol(dmat) || ncol(dmat)!=nrow(dmat))stop("Names in distance matrix do not correspond to plot names")
	row.names(dmat)=names(obj)
	names(dmat)=names(obj)
	
# find all distances from plot.tt to plot.tt + some shifter value (e.g., '3' would be plot1 to plot4, plot10 to plot3, ... plotN to plotN+3)
# tabulate these values and find the average distance from each plot to every plot+'shifter'
	torus=rep(1:ncol(obj),2)
	plus.array=array(dim=c(1,ncol(obj)))
	torus.array=array(dim=c(ncol(obj), ncol(obj)))
	for(plus in 1:ncol(obj)) {
		for(tt in 1:ncol(torus.array)) {
			from=tt
			to=torus[tt+plus]
			d.tt=dmat[from, to]
			torus.array[tt,plus]=d.tt
		}
		plus.array[1,plus]=mean(torus.array[,plus],na.rm=TRUE)
	}
	
	plus.array=as.data.frame(plus.array)
	names(plus.array)=names(obj)
	plus.array=plus.array[,order(plus.array)]
	
# randomly generate a value between 0 and maximum average distance from a plot to every other plot+'shifter'
# shift species abundances by the randomly chosen 'shifter' 
	for(ss in 1:nrow(obj)){
		if(!flag) {
			shifter=as.numeric(names(plus.array))[min(which(plus.array>=runif(1,min=min(plus.array), max=max(plus.array))))]
		} else {shifter = sample(as.numeric(names(plus.array)),1)}
		
		t.array=array(dim=c(ncol(obj),2))
		t.array[,1]=1:ncol(obj)
		for(o in 1:ncol(obj)){
			tt=torus[shifter+(o-1)]
			t.array[tt,2]=obj[ss,o]
		}
		obj[ss,]=t.array[,2]
	}
	if(flag) message("Plots were assumed to be equidistant from one another.")
	res=obj
	names(res)=names.orig
	return(res)
}

#' @title randomizing a community phylogenetics matrix: '3x' of Hardy (2008)
#'
#' @description Used for resampling data within a community dataset
#'
#' @param obj a community dataset in spacodiR format (see as.spacodi)
#' @param level numeric vector of length one giving a proportion specifying the extent of data shuffling
#' 
#' @author Jonathan Eastman
#' 
#' @details 
#' A resampling procedure for a species-by-plots matrix, based on Gotelli swapping. 
#' Shuffles abundances within species and for a pair of plots. The level defines the degree of sampling, 
#' with larger values dictating a higher level of reshuffling. For instance, if level = 0.1 and the dataset 
#' involves 20 species and 20 plots, a total of 40 (0.1x20x20) Gotelli swaps are performed.
#' 
#' @return A shuffled dataset
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' GOTELLI NJ. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606-2621.
#' 
#' @examples 
#' data(sp.example)
#' attach(sp.example)
#' spl
#' 
#' 
#' # shuffle dataset
#' resamp.3x(obj=spl, level=0.2)
#' 
#' @export
resamp.3x <-
  function(obj, level = 0.1) {
    swaps = round(level * ncol(obj) * nrow(obj))
    orig = obj
    for (swap in 1:swaps) {
      rcol = sample(1:ncol(obj), 2)
      rspp = sample(1:nrow(obj), 2)
      obj[rspp[1], rcol[1]] = orig[rspp[1], rcol[2]]
      obj[rspp[1], rcol[2]] = orig[rspp[1], rcol[1]]
      obj[rspp[2], rcol[1]] = orig[rspp[2], rcol[2]]
      obj[rspp[2], rcol[2]] = orig[rspp[2], rcol[1]]
      orig = obj
    }
    if (sum(obj) != sum(orig))
      warning("A poor result is likely.")
    return(obj)
  }


#' @title partial phylogeny randomization for tips
#'
#' @description used to shuffle tips subtended by a set of internal nodes determining by divergence time or specified by the user
#'
#' @param phy a phylogenetic tree of class phylo
#' @param node numeric value(s), specifying the internal node(s) whose tips to shuffle; see nodelabels
#' @param time.threshold numeric
#' @param proportion logical
#' 
#' @author Jonathan Eastman
#' 
#' @details
#' Either a numeric vector is supplied for node or a time.threshold. If given a set of nodes, 
#' this function will naively shuffle tips descended from the nodes in the order supplied 
#' (without regard to whether any internal node in the vector is a descendant of any other node in the node vector). 
#' If given a time.threshold, tips will be reshuffled within non-nested clades that have a rootmost node that occurs 
#' within the range [0, time.threshold]. Note that regard to absolute divergence times can be enforced with proportion=FALSE. 
#' Note further that resamp.phy(phy=phy, node=NULL, node.threshold=1, proportion=TRUE) achieves the same effect as resamp.1s.
#' 
#' @return A phylogenetic tree whose tips have been shuffled (without any modification of the underlying topology)
#' 
#' @references 
#' HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
#' 
#' @examples 
#' \dontrun{
#' data(sp.example)
#' attach(sp.example)
#' 
#' # reshuffle within a time range
#' time=1/3
#' bb <- ape::branching.times(phy)
#' bb <- bb/max(bb)
#' nodes <- (Ntip(phy)+1):max(phy$edge)
#' nodes[bb<=time] <- 1
#' dev.new()
#' plot(resamp.phy(phy, time.threshold=time, proportion=TRUE))
#' mtext("reshuffled phylogeny showing affected nodes")
#' nodelabels(cex=ifelse(nodes==1, 2, NA), col=ifelse(nodes==1, 1, NA), pch=19) 
#' }

#' 
#' @importFrom ape Ntip branching.times
#' 
#' @export
resamp.phy <-
  function(phy,
           node = NULL,
           time.threshold = 1,
           proportion = TRUE) {
    if (all(is.numeric(node))) {
      # allow tip reshuffling by prespecified node
      
      if (all(node > Ntip(phy)) & all(node <= max(phy$edge))) {
        internals = node
      } else {
        stop("Supplied 'node' value(s) appear misspecified.")
      }
    } else {
      # reshuffling of tips by temporal constraint
      
      n <- ape::Ntip(phy)
      b <- ape::branching.times(phy)
      if (proportion)
        b = b / max(b)
      tt = sapply(b, function(x)
        withinrange(x, 0, time.threshold))
      if (any(tt)) {
        tt = as.numeric(names(tt)[which(tt)])
        affected.nodes = c()
        internals = c()
        to.do = lapply(tt, function(x)
          get.descendants.of.node(x, phy, tips = FALSE))
        for (nn in 1:length(to.do)) {
          if (any(!to.do[[nn]] %in% affected.nodes)) {
            internals = c(internals, tt[nn])
          }
          affected.nodes = c(affected.nodes, to.do[[nn]][!to.do[[nn]] %in% affected.nodes])
        }
      } else {
        warning("No resampling possible")
        return(NULL)
      }
    }
    
    # SHUFFLE TIPS for NON-NESTED CLADES
    to.rand = lapply(internals, function(x)
      get.descendants.of.node(x, phy, tips = TRUE))
    for (nn in 1:length(to.rand)) {
      phy$tip.label[to.rand[[nn]]] = sample(phy$tip.label[to.rand[[nn]]])
    }
    
    return(phy)
  }


