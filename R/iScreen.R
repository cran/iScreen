#' Check WellID format
#' 
#' Examination of variable WellID format
#' @param x A character or a vector of characters as variable WellID
#' @return The value is logical of length one. 
#' Return TRUE if and only if all elements in x accord with WELLID requirement, otherwise return FALSE.
#' @examples #WellID format should be "A01", "B02", "C03", ..."J12", ...
#' WellCheck(c("A01"))
#' WellCheck(c("A1"))
#' WellCheck(c("A01", "B02", "C03"))
#' WellCheck(c("A01", "B2", "C03"))
#' @export

WellCheck <- function(x){
  x <- toupper(as.character(x))                       # converted into character vector
  if(any(nchar(x) != 3)){                             # length of character is restricted to 3
    return(FALSE)
  }else if(any(!substr(x, 1, 1) %in% LETTERS)){       # first character has to be LETTER
    return(FALSE)
  }else if(any(is.na(as.integer(substr(x, 2, 3))))){  # last two character has to be interger
    return(FALSE)
  }else{
    return(TRUE)
  }
}

#' Convert WellID
#' 
#' Convert Well ID into plot coordinate and label
#' @param WellNO A character or a vector of characters as variable WellID to converted 
#' @return A list consisting of 4 vectors. 
#' \item{x.loc}{A vector of the length of WellNO indicating x coordinate of each WellNO upon plot}
#' \item{y.loc}{A vector of the length of WellNO indicating y coordinate of each WellNO upon plot}
#' \item{x label}{A vector of names to be used as label on X axis upon plot}
#' \item{y.label}{A vector of names to be used as label on Y axis upon plot}
#' @export
WellToLoc <- function(WellNO){
  WellNO <- as.vector(toupper(as.character(WellNO)))
  if(!WellCheck(WellNO)){
    cat("Well Number Formating Wrong", "See Help For Details", sep="\n")
  }else {
    row.num <- substr(WellNO, 1, 1)                # extract row name (y axis label)
    col.num <- substr(WellNO, 2, 3)                # extract col number (x axis label)
                          
    y.label <- rev(sort(unique(row.num)))          # y axis label
    x.label <- sort(unique(col.num))               # x axis label    
                                
    y.loc <- match(row.num, LETTERS)               # y axis coordinates
    y.loc <- max(y.loc)-y.loc+1
    x.loc <- as.integer(col.num)  
    x.loc <- x.loc-min(x.loc)+1                    # x axis coordinates

    return(list(x.loc=x.loc, y.loc=y.loc, x.label=x.label, y.label=y.label))
  }
}

#' Generalized linear model for image-based HTS 
#' 
#' Perform generalized linear model for image-based HTS data analysis. 
#' @param data A data frame containing the variables in the model with the first column as Well ID. 
#' For details of the Well ID format, see \code{\link{WellCheck}}. If not found, an error will be reported.
#' @param formula An object of class \code{\link{formula}}or one that can be coerced to that class.
#' a symbolic description of the model to be fitted. 
#' The details of model specification are given under 'Details' from \code{\link{glm}}.
#' @param control A optional vector of the length as the number of rows in data indicating which row or rows
#' from data should be used as control when fitting generalized linear model. 
#' @param ... Arguments to be passed to method. 
#' For details also see \code{\link{formula}}.
#' @return A list of two compoments, fit and coefficients. 
#' \item{fit}{A fitted generlized linear model object from \code{\link{glm}}}
#' \item{coefficients}{Estimations of coefficients for each Well ID}
glm.iScreen <- function(data, formula, control=NULL, ...){
  if(is.null(data)){
    stop("argument data missing!")
  }
  if(is.null(formula)){
    stop("argument formula missing!")
  }
  if(!is.null(control)){
    if(!is.logical(control) | length(control) != nrow(data)){
      stop("argument control invalid!")
    }
  } 
  
  data$WellID <- as.character(data$WellID)
  WellID <- sort(unique(as.character(data$WellID)))
  coefficients <- data.frame(WellID=paste("WellID", WellID, sep=""))
  
  if(!is.null(control)){
    
    data.control <- data[control, ]
    data.control$WellID  <- "A00"
    data.control <- rbind(data.control, data)
    data.control$WellID <- factor(data.control$WellID)
    fit.glm <- glm(formula=formula, data=data.control, ...)
    
  }else{
    data$WellID  <- factor(data$WellID)
    fit.glm <- glm(formula=formula, data=data, ...)
  }
  
  fit.sum <- data.frame(summary(fit.glm)$coefficient)
  coefficients <- merge(coefficients, fit.sum, by.x=1, by.y=0, all.x=T, sort=F)
  coefficients$WellID <- substr(coefficients$WellID, 7, 9)
  names(coefficients)[length(names(coefficients))] <- "p.value"
  coefficients <- coefficients[order(coefficients$WellID), ]
  
  return(list(fit=fit.glm, coefficients=coefficients))
}

#' Custom analysis of image-based HTS 
#' 
#' Perform custom analysis of image-based HTS data. 
#' @param data A data frame containing the variables in the model with the first column as Well ID. 
#' For details of the Well ID format, see \code{\link{WellCheck}}. If not found, an error will be reported.
#' @param FUN A user-provided function for analysis 
#' @details For argument FUN, user has to provide a custom function. 
#' Function argument has to be from data. 
#' Return value can be a single value or a vector. 
#' For vector return, the first value has to be estimate for each Well ID. 
#' @return A list of two components, fit and coefficients. 
#' \item{fit}{NULL}
#' \item{coefficients}{Estimations of coefficients for each Well ID}
fun.iScreen <- function(data, FUN){
  if(is.null(data)){
    stop("argument data missing!")
  }
  if(is.null(FUN)){
    stop("argument FUN missing!")
  }
  
  dat.list <- split(data, data$WellID)
  res.list <- list()
  
  for(i in names(dat.list)){
    res.list[[i]] <- FUN(dat.list[[i]])
  }
  res.out <- do.call(rbind, res.list)
  coefficients <- data.frame(WellID=row.names(res.out), res.out)
  row.names(coefficients) <- NULL

  return(list(fit=NULL, coefficients=coefficients))
}

#' image-based HTS analysis
#' 
#' Analysis of image-based high-throughput RNAi screen via either generalized linear model 
#' or customized user function. 
#' @param data A data frame containing the variables in the model with the first column as Well ID. 
#' For details of the Well ID format, see \code{\link{WellCheck}}. If not found, an error will be reported.
#' @param formula An object of class "\code{\link{formula}}" (or one that can be coerced to that class): 
#' a symbolic description of the model to be fitted. 
#' The details of model specification are given under 'Details' from \code{\link{glm}}.
#' @param control A optional vector of the length as the number of rows in data indicating which row or rows
#' from data should be used as control when fitting generalized linear model. 
#' @param FUN A user-provided function for analysis 
#' @param ... Arguments to be passed to method. 
#' For details also see \code{\link{formula}}.
#' @details For argument FUN, user has to provide a custom function. 
#' Function argument has to be from data. 
#' Return value can be a single value or a vector. 
#' For vector return, the first value has to be estimate for each Well ID. 
#' @return A "iScreen" list of two components, fit and coefficients. 
#' \item{fit}{A fitted generalized linear model object from \code{\link{glm}} or NULL if FUN is used}
#' \item{coefficients}{Estimations of coefficients for each Well ID}
#' @examples
#' data(autophagy)
#' fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control  == 1))
#' head(fit.auto$coefficients)
#' @export
iScreen <- function(data=NULL, formula=NULL, control=NULL, FUN=NULL, ...){
  if(is.null(data)){
    stop("Data missing!")
  }
  if(!WellCheck(data$WellID)){
    stop("Well Number Formatting Wrong, See Help For Details")
  }
  if(is.null(formula) & is.null(FUN)){
    stop("Please specify formula or FUN, but not both")
  }
  if(!is.null(formula) & !is.null(FUN)){
    stop("Please specify formula or FUN, but not both")
  }
  if(!is.null(formula)){
    res.out <- glm.iScreen(data=data, formula=formula, control=control, ...)
  }else{
    res.out <- fun.iScreen(data=data, FUN=FUN, ...)
  }
  class(res.out) <- "iScreen"
  return(res.out)
}

#' image-based HTS Plate
#' 
#' image-based HTS Plate object to be used for analysis and visualization
#' @param data A data frame containing the variables in the model with the first column as Well ID. 
#' For details of the Well ID format, see \code{\link{WellCheck}}. If not found, an error will be reported.
#' @param column Specify the column to be used for analysis and visualization. 
#' @param log Specify if data need be logarithm transformed.
#' @param FUN A function to compute the summary statistics which can be applied to all data subsets.
#' @param ... Further arguments passed to or used by methods. See \code{\link{aggregate}} for further information.
#' @return A object "iPlate" list of 4 elements. 
#' \item{z}{A data frame containing Well ID and summary statistics for each Well ID}
#' \item{loc}{A list of coordinate information returned by \code{\link{WellToLoc}}}
#' \item{log}{A logical value indicating if data is logarithm transformed}
#' \item{FUN}{the function to compute the summary statistics}
#' @export
iPlate <- function(data, column, log=FALSE, FUN=mean, ...){

  z <- data[, column]
  
  if(log == TRUE){
    if(any(z == 0)){
      print("0 value existing and pseudocount added by 1")
      z <- aggregate(log(z+1), by=list(WellID=data$WellID), FUN=FUN, ...)
    }else{
      z <- aggregate(log(z), by=list(WellID=data$WellID), FUN=FUN, ...)
    }
  }else{
    z <- aggregate(z, by=list(WellID=data$WellID), FUN=FUN, ...)
  }
  
  loc <- WellToLoc(z$WellID)
  well.id <- data.frame(WellID=sort(paste(rep(loc$y.label, each=length(loc$x.label)), loc$x.label, sep="")))
  z <- merge(well.id, z, by.x=1, by.y=1, all.x=T, sort=F)
  z <- z[order(z$WellID), ]
  loc <- WellToLoc(z$WellID)
  
  res.out <- list(z=z, loc=loc, log=log, FUN=FUN)
  class(res.out) <- "iPlate"
  return(res.out)
  
}

#' Plotting iPlate object
#' 
#' Function for plotting object returned by iPlate. For more details about the graphical parameter arguments, see \code{\link{par}}.
#' @param object A iPlate object. 
#' @param ... Arguments to be passed to methods. See \code{\link{plot}} and \code{\link{par}}
#' @export
#' @examples
#' data(autophagy)
#' p1 <- iPlate(autophagy, "dot.number", log=TRUE, FUN=mean)
#' iPlatePlot(p1)
iPlatePlot <- function(object, ...){
  
  z <- object$z
  loc <- object$loc
  
  rgb.palette <- colorRampPalette(c("green","black","red"), space = "rgb")
  
  y <- z[, 2]
  x.label <- loc$x.label
  y.label <- loc$y.label
  cn <- length(x.label)
  rn <- length(y.label)
  
  y.len <- length(y)
  color.code.y <- y - min(y, na.rm=T)  
  color.code.y <- color.code.y/max(color.code.y, na.rm=T)  
  color.code.y <- floor(color.code.y*(y.len - 1))+1
  color.code.y <- rgb.palette(y.len)[color.code.y]  
  color.code.y <- ifelse(is.na(color.code.y), "black", color.code.y)
  
  pch.code.y <- ifelse(is.na(y), 1, 16)
  
  plot(c(0.5,cn+0.5), c(0.5,rn+0.5), type="n", 
       xaxt="n", yaxt="n", xlab="Column", ylab="Row", ...)
  points(loc$x.loc, loc$y.loc, col=color.code.y, pch=pch.code.y, ...)
  axis(1, at=1:cn, labels=x.label, ...)
  axis(2, at=1:rn, labels=y.label, ...)

}

#' Plotting iPlate legend
#' 
#' Function for plotting legend for object returned by iPlate. 
#' For more details about the graphical parameter arguments, see \code{\link{image}} and \code{\link{par}}.
#' @param object A iPlate object. 
#' @param ... Arguments to be passed to methods. See \code{\link{image}} and \code{\link{par}}.
#' @export
#' @examples
#' data(autophagy)
#' p1 <- iPlate(autophagy, "dot.number", log=TRUE, FUN=mean)
#' iPlateLegend(p1)
iPlateLegend <- function(object, ...){
  z <- object$z
  x <- z[, 2]
  ix <- 1
  nlevel=16
  minz <- min(x)
  maxz <- max(x)
  binwidth <- (maxz - minz)/nlevel
  midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  rgb.palette <- colorRampPalette(c("green","black","red"), 
                                  space = "rgb")
  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", col = rgb.palette(13), ...)
  axis(4, mgp = c(3, 1, 0), las = 2, ...)
}

#' Plotting iPlate boxplot
#' 
#' Function for box plotting object returned by iPlate. 
#' For more details about the graphical parameter arguments, see \code{\link{boxplot}} and \code{\link{par}}.
#' @param object A iPlate object. 
#' @param by Default is "row" and can also be "column".
#' @param ... Arguments to be passed to methods. See \code{\link{boxplot}} and \code{\link{par}}.
#' @export
#' @examples
#' data(autophagy)
#' p1 <- iPlate(autophagy, "dot.number", log=TRUE, FUN=mean)
#' iPlateBoxplot(p1)
iPlateBoxplot <- function(object, by="row",...){
  z <- object$z
  loc <- object$loc
  x.loc <- loc$x.loc
  y.loc <- max(loc$y.loc)-loc$y.loc+1
  x.label <- loc$x.label
  y.label <- rev(loc$y.label)
  if(by == "column"){
    boxplot(z[, 2]~x.loc, names=x.label, ...)
  }else if(by == "row"){
    boxplot(z[, 2]~y.loc, names=y.label,  ...)
  }
}

#' Plotting iScreen
#' 
#' Function for plotting object returned by iScreen. For more details about the graphical parameter arguments, see \code{\link{par}}.
#' @param object A iScreen object. 
#' @param xlab Default is "Column". 
#' @param ylab Default is "Row". 
#' @param ... Arguments to be passed to methods. See \code{\link{plot}} and \code{\link{par}}
#' @export
#' @examples
#' data(autophagy)
#' fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control == 1))
#' iScreenPlot(fit.auto)
iScreenPlot <- function(object, xlab="Column", ylab="Row", ...){
  z <- object$coefficients
  loc <- WellToLoc(z$WellID)
  well.id <- data.frame(WellID=sort(paste(rep(loc$y.label, each=length(loc$x.label)), loc$x.label, sep="")))
  z <- merge(well.id, z, by.x=1, by.y=1, all.x=T, sort=F)
  z <- z[order(z$WellID), ]
  loc <- WellToLoc(z$WellID)
  
  rgb.palette <- colorRampPalette(c("green","black","red"), space = "rgb")
  
  y <- z[, 2]
  x.label <- loc$x.label
  y.label <- loc$y.label
  cn <- length(x.label)
  rn <- length(y.label)
  
  y.len <- length(y)
  color.code.y <- y - min(y, na.rm=T)  
  color.code.y <- color.code.y/max(color.code.y, na.rm=T)  
  color.code.y <- floor(color.code.y*(y.len - 1))+1
  color.code.y <- rgb.palette(y.len)[color.code.y]  
  color.code.y <- ifelse(is.na(color.code.y), "black", color.code.y)
  
  pch.code.y <- ifelse(is.na(y), 1, 16)
  
  plot(c(0.5,cn+0.5), c(0.5,rn+0.5), type="n", 
       xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, ...)
  points(loc$x.loc, loc$y.loc, col=color.code.y, pch=pch.code.y, ...)
  axis(1, at=1:cn, labels=x.label, ...)
  axis(2, at=1:rn, labels=y.label, ...)
}

#' Plotting iScreen legend
#' 
#' Function for plotting legend for object returned by iScreen. 
#' For more details about the graphical parameter arguments, see \code{\link{image}} and \code{\link{par}}.
#' @param object A iScreen object. 
#' @param ... Arguments to be passed to methods. See \code{\link{image}} and \code{\link{par}}.
#' @export
#' @examples
#' data(autophagy)
#' fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control == 1))
#' iScreenLegend(fit.auto)
iScreenLegend <- function(object, ...){
  z <- object$coefficients
  x <- z[, 2]
  ix <- 1
  nlevel=16
  minz <- min(x)
  maxz <- max(x)
  binwidth <- (maxz - minz)/nlevel
  midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  rgb.palette <- colorRampPalette(c("green","black","red"), 
                                  space = "rgb")
  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", col = rgb.palette(13), ...)
  axis(4, mgp = c(3, 1, 0), las = 2, ...)
}

#' Plotting iScreen boxplot
#' 
#' Function for box plotting object returned by iScreen. 
#' For more details about the graphical parameter arguments, see \code{\link{boxplot}} and \code{\link{par}}.
#' @param object A iScreen object. 
#' @param by Default is "row" and can also be "column".
#' @param ... Arguments to be passed to methods. See \code{\link{boxplot}} and \code{\link{par}}.
#' @export
#' @examples
#' data(autophagy)
#' fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control == 1))
#' iScreenBoxplot(fit.auto)
iScreenBoxplot <- function(object, by="row",...){
  z <- object$coefficients
  loc <- WellToLoc(z$WellID)
  x.loc <- loc$x.loc
  y.loc <- max(loc$y.loc)-loc$y.loc+1
  x.label <- loc$x.label
  y.label <- rev(loc$y.label)
  if(by == "column"){
    boxplot(z[, 2]~x.loc, names=x.label, ...)
  }else if(by == "row"){
    boxplot(z[, 2]~y.loc, names=y.label,  ...)
  }
}

#' Batch process of image-based HTS
#' 
#' Batch processing of image-based high-throughput RNAi screen via either generalized linear model 
#' or customized user function. 
#' @param data A data frame containing the variables in the model with the first two column as Well ID and Plate ID. 
#' For details of the Well ID format, see \code{\link{WellCheck}}. If not found, an error will be reported.
#' @param formula An object of class "\code{\link{formula}}" (or one that can be coerced to that class): 
#' a symbolic description of the model to be fitted. 
#' The details of model specification are given under 'Details' from \code{\link{glm}}.
#' @param control A optional verctor of the length as the number of rows in data indicating which row or rows
#' from data should be used as control when fitting generalized linear model. 
#' @param FUN A user-provided function for anlaysis. 
#' @param ... Arguments to be passed to method. 
#' For details also see \code{\link{formula}}.
#' @details For argument FUN, user has to providde a custom function. 
#' Function argument has to be from data. 
#' Return value can be a single value or a vector. 
#' For vector return, the first value has to be estimate for each Well ID. 
#' @return A "iScreenInBatch" data frame containing summary coefficients.
#' @export
iScreenInBatch <- function(data=NULL, formula=NULL, control=NULL, FUN=NULL, ...){
  dat.list <- split(data, data$PlateID)
  res.list <- list()
  
  for(i in names(dat.list)){
    res.list[[i]] <- iScreen(dat.list[[i]], formula=formula, control=control, FUN=FUN, ...)
    res.list[[i]] <- res.list[[i]]$coefficients
    res.list[[i]] <- data.frame(PlateID=i, res.list[[i]])
  }
  
  res.out <- do.call(rbind, res.list)
  
  class(res.out) <- "iScreenInBatch"

  return(res.out)
}

#' Plotting iScreenInBatch boxplot
#' 
#' Function for plotting object returned by iScreenInBatch. For more details about the graphical parameter arguments, see \code{\link{par}}.
#' @param object A iScreenInBatch object. 
#' @param ... Arguments to be passed to methods. See \code{\link{boxplot}} and \code{\link{par}}
#' @export
iScreenInBatchBoxplot <- function(object, ...){
  boxplot(object[, 3]~object[, 1], ...)
}

#' image-based HTS Well object
#' 
#' image-based HTS Well object to be used for visualization
#' @param x A vector containing X coordinate of each plot unit.
#' @param y A vector containing Y coordinate of each plot unit.
#' @param d A vector containing diameter of each plot unit. 
#' @param c A vector containing plot colore of each plot unit.
#' @param n A numeric or vector containing number of sides for each plot unit.
#' @param angle A numeric or vector containing rotation angle of each plot unit, in degrees.
#' @param type A numeric or vector containing plotting type type=1 => interior filled, type=2 => edge, type=3 => both.
#' @return A object "iWell". 
#' @details For details also see \code{\link{ngon}}. 
#' @export
iWell <- function(x, y, d, c, n=4, angle=0, type=1){
  stopifnot(is.numeric(x), is.numeric(y), is.numeric(d), is.numeric(n), 
            is.numeric(angle), type %in% c(1, 2, 3))
  if(is.factor(c)){
    c <- as.character(c)
  }
  dat <- list(df=data.frame(x=x, y=y, d=d, c=c), n=n, angle=angle, type=type)
  class(dat) <- "iWell"
  return(dat)
}

#' Plotting iWell
#' 
#' Function for plotting object returned by iPlate. For more details about the graphical parameter arguments, see \code{\link{par}}.
#' @param object A iScreen object. 
#' @param xlab Default is "X". 
#' @param ylab Default is "Y". 
#' @param ... Arguments to be passed to methods. See \code{\link{plot}} and \code{\link{par}}
#' @import maptree
#' @export
iWellPlot <- function(object, xlab="X", ylab="Y", ...){
  
  plot(object$df$x, object$df$y, type="n", xlab=xlab, ylab=ylab, ...)
  apply(object$df, 1, function(x){
    ngon(x, object$n, object$angle, object$type)
  })
}
