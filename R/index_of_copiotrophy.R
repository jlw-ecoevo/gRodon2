
#' Create an Index of Copiotrophy
#'
#' This function creates an index of copiotrophy from a set of genomes that captures
#' the range of growth strategies encoded in that set.
#' (as in https://doi.org/10.1101/2025.09.01.673550 and https://doi.org/10.1126/science.ado5323).
#'
#' @param x Dataframe with the columns "dCUB", "nCAZy", and "nGenes" containing dCUB
#' values from gRodon's predictGrowth() function, the number of CAZymes per genome
#' (as predicted by dbCAN or other tool), and the number of protein coding genes per genome
#' respectively. If x is omitted, these values can be supplied as vectors.
#' @param dCUB Vector of dCUB values from gRodon's predictGrowth() function. Only required
#' if x is not provided as an argument.
#' @param nCAZy Vector of the number of CAZymes per genome. Only required
#' if x is not provided as an argument.
#' @param nGenes Vector of the number of protein coding genes per genome. Only required
#' if x is not provided as an argument.
#' @return fitIoC returns a list with the following elements:
#' \describe{
#'   \item{IoC}{The calculated IoC values for each of the genomes in x}
#'   \item{model}{The fitted IoC model that can be passed to the predictIoC() function
#'   to place new genomes on this fitted index of copiotrophy.}
#' }
#' @examples
#' #First let's take a look at the "mags" dataframe provided in this package
#' # (from https://doi.org/10.1101/2025.09.01.673550)
#' load(system.file('extdata',
#'   'mags_df.rda',
#'   package = 'gRodon'))
#' head(mags)
#'
#' #Let's hold out 3 mags for illustration
#' mags_held_out <- mags[1:3,]
#' mags_train <- mags[4:nrow(mags),]
#'
#' #Now let's fit an index of copiotrophy using these mags
#' IoC <- fitIoC(mags_train)
#'
#' #Now let's place our held out mags onto this fitted IoC
#' predictIoC(mags_held_out,model=IoC)
#'
#' @export
fitIoC <- function(x=NULL,dCUB=NA,nCAZy=NA,nGenes=NA){
  if(is.null(x) | class(x)!="data.frame"){
    if(is.na(dCUB+nCAZy+nGenes)){
      stop("If not providing a dataframe \"x\" with the columns, \"dCUB\",\"nCAZy\", and \"nGenes\", then please provide these values as vector arguments")
    }
    else if(length(dCUB)!=length(nCAZy) | length(nCAZy)!=length(nGenes)){
      stop("The arguments \"dCUB\",\"nCAZy\", and \"nGenes\" must have the same length")
    }

    x <- data.frame(dCUB=as.numeric(dCUB),
                    nCAZy=as.numeric(nCAZy),
                    nGenes=as.numeric(nGenes))

  } else if(sum(c("dCUB","nCAZy","nGenes") %in% names(x))<3){
    stop("\"x\" must have the columns, \"dCUB\",\"nCAZy\", and \"nGenes\"")
  }

  x$rCAZy <- x$nCAZy/x$nGenes
  x.pca <- prcomp(x[,c("rCAZy","dCUB","nGenes")],scale=T)
  x$PC1 <- x.pca$x[,1]
  x$IoC <- x$PC1-min(x$PC1)+1

  return(list(IoC=x$IoC,model=x.pca))
}


#' Calculate Index of Copiotrophy for a Set of Genomes
#'
#' This function calculates IoC values, given some already fitted index of copiotrophy
#' (as in https://doi.org/10.1101/2025.09.01.673550 and https://doi.org/10.1126/science.ado5323).
#'
#' @param x Dataframe with the columns "dCUB", "nCAZy", and "nGenes" containing dCUB
#' values from gRodon's predictGrowth() function, the number of CAZymes per genome
#' (as predicted by dbCAN or other tool), and the number of protein coding genes per genome
#' respectively. If x is omitted, these values can be supplied as vectors.
#' @param dCUB Vector of dCUB values from gRodon's predictGrowth() function. Only required
#' if x is not provided as an argument.
#' @param nCAZy Vector of the number of CAZymes per genome. Only required
#' if x is not provided as an argument.
#' @param nGenes Vector of the number of protein coding genes per genome. Only required
#' if x is not provided as an argument.
#' @param model The output of fitIoC(), which creates an index of copiotrophy from a set of genomes that captures
#' the range of growth strategies encoded in that set. Alternatively, set to "permafrost" to use the IoC
#' computed in https://doi.org/10.1101/2025.09.01.673550 from MAGs from a permafrost warming experiment
#' or "pacific" to use the IoC from the P16 Pacific ocean transect in https://doi.org/10.1126/science.ado5323
#' @return predictIoC returns a vector of IoC values
#' @examples
#' #First let's take a look at the "mags" dataframe provided in this package
#' # (from https://doi.org/10.1101/2025.09.01.673550)
#' load(system.file('extdata',
#'   'mags_df.rda',
#'   package = 'gRodon'))
#' head(mags)
#'
#' #Let's hold out 3 mags for illustration
#' mags_held_out <- mags[1:3,]
#' mags_train <- mags[4:nrow(mags),]
#'
#' #Now let's fit an index of copiotrophy using these mags
#' IoC <- fitIoC(mags_train)
#'
#' #Now let's place our held out mags onto this fitted IoC
#' predictIoC(mags_held_out,model=IoC)
#'
#' #Alternatively, use one of the pre-computed IoC models available in this package
#' predictIoC(mags_held_out,model="permafrost")
#' predictIoC(mags_held_out,model="pacific")
#'
#' @export
predictIoC <- function(x=NULL,dCUB=NA,nCAZy=NA,nGenes=NA,model){
  if(class(model)=="list"){
    if(class(model$model)!="prcomp"){
      stop("Please provide a valid IoC model either from the fitIoC() or prcomp() functions, or one of the built-in models: \"permafrost\" or \"pacific\"")
    } else {
      model <- model$model
    }
  }

  if(class(model)=="character"){
    if(model=="permafrost"){
      model <- IoC.permafrost
    } else if(model=="pacific") {
      model <- IoC.pacific
    }
  } else if(class(model)!="prcomp"){
    stop("Please provide a valid IoC model either from the fitIoC() or prcomp() functions, or one of the built-in models: \"permafrost\" or \"pacific\"")
  }

  if(is.null(x) | class(x)!="data.frame"){
    if(is.na(dCUB+nCAZy+nGenes)){
      stop("If not providing a dataframe \"x\" with the columns, \"dCUB\",\"nCAZy\", and \"nGenes\", then please provide these values as vector arguments")
    }
    else if(length(dCUB)!=length(nCAZy) | length(nCAZy)!=length(nGenes)){
      stop("The arguments \"dCUB\",\"nCAZy\", and \"nGenes\" must have the same length")
    }

    x <- data.frame(dCUB=as.numeric(dCUB),
                    nCAZy=as.numeric(nCAZy),
                    nGenes=as.numeric(nGenes))

  } else if(sum(c("dCUB","nCAZy","nGenes") %in% names(x))<3){
    stop("\"x\" must have the columns, \"dCUB\",\"nCAZy\", and \"nGenes\"")
  }

  x$rCAZy <- x$nCAZy/x$nGenes
  x.pred <- predict(model,x[,c("rCAZy","dCUB","nGenes")])
  x$PC1 <- x.pred[,1]
  x$IoC <- x$IoC <- x$PC1-min(model$x[,1])+1

  return(x$IoC)
}


