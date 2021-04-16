#' @title Internal function to check matrix structure
#'
#' @description Check matrix structure
#'
#' @param obj object to check
#' @param tol phylo object
#' @author JM Eastman
#' 
#' @return a matrix of phylogenetic distance for all pairs of species
#' 
#' 
check.distmat <- function(obj, tol = 1e-9) {
  if (any(grepl("matrix", class(obj)))) {
    if (all(diag(obj) - 0 <= tol)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}

#' @title Internal function find undersampled plots
#'
#' @description find undersampled plots
#'
#' @param sp.plot species sample matrix
#' @param verbose logical
#' 
#' @author JM Eastman
#' 
#' @return a species sample matrix
#' 
prune.sp = function(sp.plot, verbose = FALSE) {
  rr = rownames(sp.plot)
  l.spp = nrow(sp.plot)
  drop.plots = vector()
  for (sp in 1:ncol(sp.plot)) {
    l.nulls = length(which(sp.plot[, sp] == 0))
    if ((l.spp - l.nulls) < 2) {
      drop.plots = cbind(drop.plots, sp)
    }
  }
  
  plot.names.orig = colnames(sp.plot)
  dropped.plots = plot.names.orig[drop.plots]
  if (length(drop.plots) != 0) {
    if (verbose)
      message({
        cat("\nThe following plots were dropped from sp.plot:\n\t")
        cat(dropped.plots, sep = " ")
        cat("\n")
      })
    sp.plot = as.matrix(sp.plot[, -as.numeric(drop.plots[!is.na(drop.plots)])])
    rownames(sp.plot) = rr
    colnames(sp.plot) = plot.names.orig[which(!plot.names.orig %in% plot.names.orig[drop.plots])]
    if (is.null(ncol(sp.plot)))
      sp.plot = NULL
  }
  return(sp.plot)
}

#' @title Internal function that Compute phylogenetic distances among species
#'
#' @description  Utility for automating the process of pruning datasets for use in spacodiR
#'
#' @param sp.plot species sample matrix
#' @param phy a phylogenetic tree of class phylo
#' @param sp.traits species traits matrix
#' @param prune logical
#' @param verbose logical
#' 
#' @author JM Eastman
#' 
#' @details 
#' Note that nearly all spacodiR functions require that trees have present (at least) all species sampled in the community dataset
#' If prune=TRUE, the list of plots pruned will be printed to the console
#' 
#' @return A list of pruned dataset(s)
#' \itemize{
#'   \item sp.plot
#'   \item sp.tree
#'   \item sp.traits
#' }
#' 
#' @examples 
#' 
#' # load a species-by-plots matrix, along with a tree
#' data(sp.example)
#' attributes(sp.example)
#' attach(sp.example)
#' spl
#' phy
#' trt
#' 
#' # prune out undersampled plots
#' spl[,2] <- 0
#' match.spacodi.data(spl) -> sp.plot.new
#' sp.plot.new
#' 
#' # match datasets where sp.traits is smaller than the remainder
#' match.spacodi.data(sp.plot = spl, phy = phy, sp.traits = trt[1:6,])
#' 
#' 
#' @export
match.spacodi.data <-
  function(sp.plot,
           phy = NULL,
           sp.traits = NULL,
           prune = TRUE,
           verbose = FALSE) {
    # major error checking
    sp.plot <- as.matrix(sp.plot)
    if (is.null(row.names(sp.plot)))
      stop("Check that sp.plot has row names as species.")
    if (is.null(colnames(sp.plot))) {
      warning("sp.plot does not appear to have plots as column names.")
      names(sp.plot) <- paste("plot", 1:ncol(sp.plot), sep = "")
    }
    
    if (!missing(phy)) {
      if (class(phy) == "phylo") {
        if (length(unique(phy$tip.label)) != length(phy$tip.label))
          stop("Redundant taxa were found in tree.")
        if (!is.null(phy$node.label))
          phy$node.label <- NULL
      }
    }
    
    if (length(unique(names(sp.plot))) != length(names(sp.plot)))
      stop("Redundant plots were found in sp.plot.")
    
    # check poor values in sp.plot
    if (any(!is.finite(sp.plot))) {
      poor <- which(!is.finite(sp.plot))
      poor %% nrow(sp.plot) -> poor.rows
      ceiling(poor / nrow(sp.plot)) -> poor.cols
      poor.data = as.list(paste(rownames(sp.plot)[poor.rows], colnames(sp.plot)[poor.cols]))
      stop(
        "Poor data values found in supplied sp.plot (NA, NaN, or Inf):\n\n\t",
        toString(poor.data),
        "\n"
      )
    }
    
    
    if (prune)
      sp.plot <- prune.sp(sp.plot, verbose = verbose)
    
    # match ordering and size of all data objects
    usable.species <- rownames(sp.plot)
    
    # phy only
    if (missing(sp.traits) & !missing(phy)) {
      if (class(phy) == "phylo") {
        usable.species.tmp <- intersect(phy$tip.label, usable.species)
        phy <- reorderspacodiobj(phy, usable.species.tmp)
        usable.species <- usable.species.tmp[match(phy$tip.label, usable.species.tmp)]
        if (prune)
          sp.plot <- prune.sp(reorderspacodiobj(obj = sp.plot, names = usable.species), verbose =
                               verbose)
        else
          sp.plot <- reorderspacodiobj(sp.plot, usable.species)
        r.out <- list(sp.plot = sp.plot, sp.tree = phy)
      } else if (check.distmat(phy)) {
        usable.species.tmp <- intersect(rownames(phy), usable.species)
        phy <- reorderspacodiobj(phy, usable.species.tmp)
        usable.species <- usable.species.tmp[match(rownames(phy), usable.species.tmp)]
        if (prune)
          sp.plot <- prune.sp(reorderspacodiobj(sp.plot, usable.species), verbose =
                               verbose)
        else
          sp.plot <- reorderspacodiobj(sp.plot, usable.species)
        r.out <- list(sp.plot = sp.plot, sp.tree = phy)
      }
    }
    
    # sp.traits only
    if (!missing(sp.traits) & missing(phy)) {
      usable.species <- intersect(usable.species, rownames(sp.traits))
      if (prune)
        sp.plot <- prune.sp(reorderspacodiobj(sp.plot, usable.species), verbose =
                             verbose)
      else
        sp.plot <- reorderspacodiobj(sp.plot, usable.species)
      r.out <- list(sp.plot = sp.plot,
                   sp.traits = reorderspacodiobj(sp.traits, usable.species))
    }
    
    # both sp.traits and phy
    if (!missing(sp.traits) & !missing(phy)) {
      if (class(phy) == "phylo") {
        usable.species.tmp <- intersect(phy$tip.label, intersect(rownames(sp.traits), usable.species))
        phy <- reorderspacodiobj(phy, usable.species.tmp)
        usable.species <- usable.species.tmp[match(phy$tip.label, usable.species.tmp)]
        if (prune) {
          sp.plot <- prune.sp(reorderspacodiobj(sp.plot, usable.species), verbose =
                                verbose)          
        } else {
          sp.plot <- reorderspacodiobj(sp.plot, usable.species)          
        }
        r.out <- list(
          sp.plot = sp.plot,
          sp.tree = reorderspacodiobj(phy, usable.species),
          sp.traits = reorderspacodiobj(sp.traits, usable.species)
        )
      } else {
        if (check.distmat(phy)) {
          usable.species.tmp <- intersect(rownames(phy), intersect(rownames(sp.traits), usable.species))
          phy <- reorderspacodiobj(phy, usable.species.tmp)
          usable.species <- usable.species.tmp[match(rownames(phy), usable.species.tmp)]
          if (prune) {
            sp.plot <- prune.sp(reorderspacodiobj(sp.plot, usable.species), verbose =
                                  verbose)            
          } else {
            sp.plot <- reorderspacodiobj(sp.plot, usable.species)            
          }
          
          r.out <- list(
            sp.plot = sp.plot,
            sp.tree = reorderspacodiobj(phy, usable.species),
            sp.traits = reorderspacodiobj(sp.traits, usable.species)
          )
        }
      } 
    }
    
    # neither sp.traits nor phy
    if (missing(sp.traits) & missing(phy)) {
      if (prune)
        sp.plot <- prune.sp(reorderspacodiobj(sp.plot, usable.species), verbose =
                             verbose)
      else
        sp.plot <- reorderspacodiobj(sp.plot, usable.species)
      r.out <- list(sp.plot = sp.plot)
    }
    
    return(r.out)
  }


#' @title Internal function find re ordering dataset
#'
#' @description re ordering dataset
#'
#' @param obj phylo object
#' @param names character vector of species names
#' 
#' @author JM Eastman
#' 
#' @importFrom ape drop.tip
#' 
#' @return a phylo object
reorderspacodiobj = function(obj, names) {
  if (any(grepl("phylo", class(obj)))) {
    obj.labels = obj$tip.label
    if (any(!obj.labels %in% names))
      obj <- ape::drop.tip(obj, obj.labels[!obj.labels %in% names])
    else
      obj <- obj
  } else if (check.distmat(obj)) {
    obj.labels <- rownames(obj)
    obj <- as.matrix(obj[match(names, obj.labels), match(names, obj.labels)])
    rownames(obj) <- colnames(obj) <- names
  } else {
    nn <- colnames(obj)
    obj.labels <- rownames(obj)
    obj <- as.data.frame(obj[match(names, obj.labels), ])
    names(obj) <- nn
  }
  if (!all(names %in% obj.labels))
    warning(paste(
      paste(names[!names %in% obj.labels], sep = " ", collapse = " "),
      "were not found in the supplied object",
      sep = " "
    ))
  return(obj)
}


#' @title Generates community phylogenetic datasets
#'
#' @description Generates community phylogenetic datasets to be used in the external program SPACoDi
#'
#' @param sp.plot a community dataset formatted for the R-package spacodiR
#' @param outfile a formatted file for the Windows executable SPACoDi
#' 
#' @details This utility writes a species-by-plots matrix into a format readable by the external program SPACoDi, 
#' a Windows executable
#' 
#' @examples 
#' # generate a community-phylogenetics dataset
#' data(sp.example)
#' attach(sp.example)
#' 
# save the dataset to working directory
#' write.spacodi.data(sp.plot=spl, outfile="spacodi.formatted.txt")
#' 
#' @importFrom utils write.table
#' 
#' @author JM Eastman
#' 
#' @return a phylo object
#' @export
write.spacodi.data <-
  function(sp.plot, outfile) {
    if (file.exists(outfile)) {
      warning("Overwrote existing outfile.")
      unlink(outfile)
    }
    names = names(sp.plot)
    for (n in 1:length(names)) {
      cat(c("\t", names[n]),
          file = outfile,
          append = TRUE,
          sep = "")
    }
    cat("\n",
        file = outfile,
        append = TRUE,
        sep = "")
    utils::write.table(
      sp.plot,
      outfile,
      quote = FALSE,
      col = FALSE,
      append = TRUE,
      sep = "\t"
    )
  }


#' @title Converts from spacodi or picante data formats into phylocom format
#'
#' @description Converts from spacodi or picante data formats into phylocom format
#'
#' @param data species-by-plots matrix
#' @param picante logical
#' @param outfile an optional text file to which to write output
#' 
#' @details This utility converts a species-by-plots matrix into triplet format, 
#' which is readable by the external program phylocom. If picante=TRUE, 
#' the data are expected to be in the form used for picante (i.e., a plots-by-species matrix; picante-package).
#'  If the user selects picante=FALSE, the data are expected to be in the form used for spacodiR (i.e., a species-by-plots matrix). 
#'  The user has the option to save an output file, defined by outfile.
#'  
#' @references 
#'  WEBB CO, DD ACKERLY and SW KEMBEL. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 24:2098-2100.
#' 
#' @examples 
#' # call example data from SPACoDi
#' data(sp.example)
#' attach(sp.example)
#' spl->d.spacodi  
#' d.spacodi ## SPACoDi format
#' 
#' # convert to phylocom
#' as.phylocom(data=spl, picante=FALSE)->d.phylocom
#' d.phylocom ## phylocom format
#' 
#' # convert dataset to picante
#' as.picante(data=d.phylocom)->d.picante
#' d.picante ## picante format
#' 
#' # convert back to SPACoDi 
#' as.spacodi(data=d.picante)
#' 
#' 
#' @author JM Eastman
#' 
#' @return A named array, formatted for use in phylocom; note that while the R-object returned by this function has column names, 
#' if output is written to a file, 
#' the header is dropped (as appropriate for use in the external phylocom executable: http://www.phylodiversity.net/phylocom/).
#' 
#' @export
as.phylocom <- function(data,
                        picante = FALSE,
                        outfile = NULL) {
  if (picante)
    data = as.spacodi(data)
  
  if (class(data) != "data.frame") {
    if (class(data) == "matrix")
      data = as.data.frame(data)
    else if (class(data) == "vector")
      data = as.data.frame(data)
    else
      stop("Function cannot handle data that are not in data.frame format")
  }
  if (length(unique(row.names(data))) != nrow(data))
    warning("Data do not appear to be in proper format.  Results may be nonsensical.")
  plots = ncol(data)
  spp = nrow(data)
  out = array(dim = c(spp * plots, 3))
  for (plot in 1:plots) {
    start = ((plot - 1) * spp) + 1
    end = plot * spp
    out[start:end, ] = cbind(names(data)[plot], data[, plot], row.names(data))
  }
  out = as.data.frame(out)
  names(out) = c("plot", "samples", "species")
  if (!is.null(outfile)) {
    utils::write.table(
      out,
      outfile,
      quote = FALSE,
      row = FALSE,
      col = FALSE,
      sep = "\t"
    )
  }
  return(out)
}


#' @title Converts from spacodi or phylocom data formats into picante format
#'
#' @description Converts from spacodi or phylocom data formats into picante format
#'
#' @param data a community dataset in phylocom or spacodi format
#' @param outfile an optional text file to which to write output
#' 
#' @details This utility converts a community dataset (either from phylocom or spacodiR) 
#' into a format interpretable by picante (see picante-package). phylocom format is also 
#' referred to as triplet-formatting, where plots are within the first column, abundances 
#' in the second, and species names in the third column of the dataframe. The user has the 
#' option to save an output file, defined by outfile. SPACoDi format is similar to that for picante, 
#' where dataframes between these packages are transposed. SPACoDi format should have species as row names.
#' 
#' 
#' @examples 
#' # call example data from SPACoDi
#' data(sp.example)
#' attach(sp.example)
#' spl->d.spacodi  
#' 
#' 
#' 
#' @author JM Eastman
#' 
#' @return An array, formatted for use in picante
#' 
#' @export
as.picante <-
  function(data, outfile = NULL) {
    if (ncol(data) == 3) {
      message("Formatting assumed to be that used with phylocom")
      if (class(data) != "data.frame") {
        if (class(data) == "matrix")
          data = as.data.frame(data)
        else if (class(data) == "vector")
          data = as.data.frame(data)
        else
          stop("Function cannot handle data that are not in data.frame format")
      }
      species = as.character(unique(data[, 3]))
      dd = split(data, data[, 1])
      out = array(dim = c(length(species) -> spp, length(dd) -> plots))
      for (plot in 1:plots) {
        cur.array = array(dim = c(spp, 1))
        ind = as.numeric(as.vector(dd[[plot]][, 2]))
        dd[[plot]] -> cur.plot
        names(ind) = as.character(cur.plot[, 3])
        for (r in 1:nrow(cur.plot)) {
          cur.array[which(names(ind[r]) == species)] = ind[r]
        }
        cur.array[which(is.na(cur.array))] = 0
        out[, plot] = as.numeric(cur.array)
      }
      out = as.data.frame(t(out))
      names(out) = species
      row.names(out) = names(dd)
      if (!is.null(outfile)) {
        write.spacodi.data(out, outfile)
      }
      return(out)
    } else {
      message("Formatting assumed to be that used with spacodi")
      if (class(data) != "data.frame") {
        if (class(data) == "matrix")
          data = as.data.frame(data)
        else if (class(data) == "vector")
          data = as.data.frame(data)
        else
          stop("Function cannot handle data that are not in data.frame format")
      }
      if (is.null(row.names(data)))
        row.names(data) = paste("plot", seq(1:nrow(data)))
      if (is.null(names(data)))
        names(data) = paste("species", seq(1:nrow(data)))
      
      out = as.matrix(t(data))
      if (!is.null(outfile)) {
        utils::write.table(out, outfile, quote = FALSE)
      }
      return(out)
    }
  }


#' @title Cnverts from picante or phylocom data formats into spacodi format.
#'
#' @description Converts from picante or phylocom data formats into spacodi format.
#'
#' @param data a community dataset in phylocom or spacodi format
#' @param outfile an optional text file to which to write output
#' 
#' @details This utility converts a community dataset (either from phylocom or picante (see picante) into a 
#' format interpretable by either this R-package spacodiR or the external program SPACoDi, a Windows executable
#' (available at http://ebe.ulb.ac.be/ebe/Software.html). Note also that the community-dataset format used here 
#' is also that called for by the package vegan; see vegandocs
#' 
#' phylocom format is also referred to as triplet-formatting, where plots are within the first column, 
#' abundances in the second, and species names in the third column of the dataframe. picante format is 
#' simply the transpose of spacodiR-formatting of the community dataset: in spacodiR, species are expected 
#' as the row names of the dataframe, where plots are represented as the column names. The user has the option 
#' to save an output file, defined by outfile.
#' 
#' @examples 
#' # call example data from spacodiR
#' data(sp.example)
#' attach(sp.example)
#' spl->d.spacodi  
#' d.spacodi ## SPACoDi format
#' 
#' @author JM Eastman
#' 
#' @return An array, formatted for use in spacodiR
#' 
#' @export
as.spacodi <- function(data, outfile = NULL) {
  if (ncol(data) != 3) {
    message("Formatting assumed to be that used with picante")
    if (class(data) != "data.frame") {
      if (class(data) == "matrix")
        data = as.data.frame(data)
      else if (class(data) == "vector")
        data = as.data.frame(data)
      else
        stop("Function cannot handle data that are not in data.frame format")
    }
    
    if (is.null(row.names(data)))
      row.names(data) = paste("plot", seq(1:nrow(data)))
    if (is.null(names(data)))
      names(data) = paste("species", seq(1:nrow(data)))
    
    out = as.matrix(t(data))
    if (!is.null(outfile)) {
      utils::write.table(out, outfile, quote = FALSE)
    }
    return(out)
  } else {
    message("Formatting assumed to be that used with phylocom")
    if (class(data) != "data.frame") {
      if (class(data) == "matrix")
        data = as.data.frame(data)
      else if (class(data) == "vector")
        data = as.data.frame(data)
      else
        stop("Function cannot handle data that are not in data.frame format")
    }
    species = as.character(unique(data[, 3]))
    dd = split(data, data[, 1])
    out = array(dim = c(length(species) -> spp, length(dd) -> plots))
    for (plot in 1:plots) {
      cur.array = array(dim = c(spp, 1))
      ind = as.numeric(as.vector(dd[[plot]][, 2]))
      dd[[plot]] -> cur.plot
      names(ind) = as.character(cur.plot[, 3])
      for (r in 1:nrow(cur.plot)) {
        cur.array[which(names(ind[r]) == species)] = ind[r]
      }
      cur.array[which(is.na(cur.array))] = 0
      out[, plot] = as.numeric(cur.array)
    }
    out = as.data.frame(out)
    names(out) = names(dd)
    row.names(out) = species
    if (!is.null(outfile)) {
      write.spacodi.data(out, outfile)
    }
    return(out)
  }
}
