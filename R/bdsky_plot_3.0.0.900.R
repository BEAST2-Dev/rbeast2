#' Plot Birth_Death Skyline Plots
#'
#' This function plots the epidemiological parameters of a BDSKY analysis with
#' their HPD intervals.
#'
#' @param logs (character) The path to the log file(s) or the path to the file with
#'              ".txt" extension that stores the path to all log files to be
#'              analyzed or NULL to select (default: NULL)
#' @param burninpercent (interger) The percentage of samples that should be ignored
#'                      from the log file
#' @param recent (numeric) Date of the most recent sample (for plotting from past
#'                to present)
#' @param gridSize (interger)
#' @param RepNumb (character) Name of the R0 parameter (default: "R0.s")
#' @param bUninfectiousRate (character) Name of the bacomeUninfectious parameter
#'           (default: "becomeUninfectiousRate.s")
#' @param sProportion (character) Name of the samplingProportion parameter
#'           (default: "samplingProportion.s")
#' @param startSampling (numeric) If samplingProportion was fixed to zero before
#'           a sampling date and >0 afterwards, the date when sampling started
#'           (default: NULL)
#' @author Denise Kuehnert (denise.kuehnert@gmail.com)
#' @author Carlo Pacioni (carlo.pacioni@gmail.com)
#' @return A PDF with the plots
#' @import s20x
#' @import boa
#' @import Hmisc
#' @import miscTools
#' @importFrom grDevices dev.copy2pdf
#' @importFrom graphics abline lines plot polygon
#' @importFrom stats median
#' @importFrom utils read.table tail
#' @export

bdsky_plot <- function(logs=NULL, burninpercent=10, recent=NULL, gridSize=20,
                       RepNumb="R0.s",
                       bUninfectiousRate="becomeUninfectiousRate.s",
                       sProportion="samplingProportion.s",
                       startSampling=NULL) {

  #### Function helper ####
  # input : a matrix M and a column ascii name
  # output : the numeric index of the column
  colnameindex <- function(M, colname0) {
    colsnames=names(M[1,]);
    theindex=which(colsnames==colname0);
    return(theindex);
  }
  # Return the extension from a file name
  getextension <- function(x) {
    substr(x, nchar(x) - 3, nchar(x))
  }
  # Read a text file and returns the loglist
  getLogNames <- function(logs) {
    loglist <- read.table(logs, as.is=TRUE, header=FALSE)
    dir.in <- dirname(logs)
    loglist[, 1] <- paste(dir.in, loglist[, 1], sep="/")
    return(loglist)
  }

  # Return a df with on character column of logfile names from a character vector
  dfLogNames <- function(logs) {
    loglist <- data.frame(logs)
    loglist[, 1] <- as.character(loglist[, 1])
    return(loglist)
  }

  #########################


  if(is.null(logs)) {
    message("Please, select the file that stores the path to all log files to be
            analyzed (with extension .txt) or the actual log file.")
    logs <- file.choose()
  }

  if(length(logs) == 1) {
    if(!file.exists(logs))  stop("Something went wrong.
              Please check that logs is a character vector with either the path to
              the file with .txt extension that stores the path to all log files
              to be analyzed, the path to the log file(s) or NULL to select the
              file interactively")
    if(getextension(logs) == ".txt") loglist <- getLogNames(logs)
    if(getextension(logs) == ".log") loglist <- dfLogNames(logs)
    } else {
      loglist <- dfLogNames(logs)
    }

  if(is.null(recent)) stop("Please provide an integer value for recent")

  burninpercent <- as.integer(burninpercent)
  gridSize <- as.integer(gridSize)

  for(i in 1:length(loglist[,1])){

    # /* read and assign file from log list */
    assign(paste("log", i, sep=''), read.table(loglist[i,], header=T))
    attach(get(paste("log", i, sep='')))

    R0_names <- names(get(paste("log", i, sep='')))[which(
               regexpr(RepNumb, names(get(paste("log", i, sep=''))))>0)]

    delta_names <- names(get(paste("log", i, sep='')))[which(
      regexpr(bUninfectiousRate, names(get(paste("log", i, sep=''))))>0)]

    sampling_names <- names(get(paste("log", i, sep='')))[which(
      regexpr(sProportion, names(get(paste("log", i, sep=''))))>0)]

    treeheights <- get(paste("log", i, sep=''))[, match(
                   "treeheight", tolower(names(get(paste("log", i, sep='')))))]
    origins <- get(paste("log", i, sep=''))$origin

    nsamples <- length(get(R0_names[1]))
    burnin <- round(burninpercent*nsamples/100)
    width <- median(origins[burnin:nsamples])

    F_intervalNumber <- length(R0_names)
    G_intervalNumber <- length(delta_names)
    H_intervalNumber <- length(sampling_names)

    if (max(F_intervalNumber, G_intervalNumber, H_intervalNumber) > gridSize)
    {gridSize <- max(F_intervalNumber, G_intervalNumber, H_intervalNumber)}

    medians <- matrix(data=NA, nrow=1, ncol=gridSize)
    medians_G <- matrix(data=NA, nrow=1, ncol=gridSize)
    medians_H <- matrix(data=NA, nrow=1, ncol=gridSize)

    hpd_F <- matrix(data=NA, nrow=2, ncol=gridSize)
    hpd_G <- matrix(data=NA, nrow=2, ncol=gridSize)
    hpd_H <- matrix(data=NA, nrow=2, ncol=gridSize)

    F <- matrix(data=NA, nrow=nsamples-burnin, ncol=gridSize) #R0
    G <- matrix(data=NA, nrow=nsamples-burnin, ncol=gridSize) #becomeuninfectiousRate
    H <- matrix(data=NA, nrow=nsamples-burnin, ncol=gridSize) #samplingProportion

    step <- width/(gridSize - 1)
    F_times <- seq(recent - width, recent, step)

    for(k in 1:(nsamples - burnin)){

      time <- origins[k + burnin]

      for (l in 1:length(F_times)){
        F_index <- ceiling(F_intervalNumber -
                            (recent - F_times[l])/(time/F_intervalNumber))
        G_index <- ceiling(G_intervalNumber -
                            (recent - F_times[l])/(time/G_intervalNumber))
        if(is.null(startSampling)) {
          H_index <- ceiling(H_intervalNumber -
                               (recent - F_times[l])/(time/H_intervalNumber))
        } else {
          if(F_times[l] > startSampling) {
            H_index <- 2
          } else {
            H_index <- 1
          }
        }


        F[k,l] <- get(R0_names[max(F_index, 1)])[k + burnin]
        G[k,l] <- get(delta_names[max(G_index, 1)])[k + burnin]
        H[k,l] <- get(sampling_names[max(H_index, 1)])[k + burnin]

      }
    }

    for(j in 1:gridSize) {
      if (length(which(F[,j]!="NA")) > (nsamples/10)) {
        medians[1,j] <- median(F[,j],na.rm=T)
        medians_G[1,j] <- median(G[,j],na.rm=T)
        medians_H[1,j] <- median(H[,j],na.rm=T)
        hpd_F[,j] <- boa.hpd(F[which(F[,j]!="NA"),j], 0.05)[1:2]
        hpd_G[,j] <- boa.hpd(G[which(G[,j]!="NA"),j], 0.05)[1:2]
        hpd_H[,j] <- boa.hpd(H[which(H[,j]!="NA"),j], 0.05)[1:2]
      }
    }

    layout20x(3,1)

    # /* plot R0 */
    plot(1, ylab=expression(R[0]), xlim=c(recent-width, recent),
         ylim=c(0,max(hpd_F[2,],na.rm=T)*1.1), xlab="Year", col='white', main='')
    minor.tick(nx=5, ny=2, tick.ratio=.2)
    polygon(c(F_times,rev(F_times)), c(hpd_F[2,], rev(hpd_F[1,])), col="grey90",
            border = NA)
    lines(c(F_times), c(medians[1,]), type='l')
    abline(1,0,col='grey')

    # /* plot become non-infectious rate */
    plot(1, ylab=expression(delta), xlim=c(recent-width, recent),
         ylim=c(0,max(hpd_G[2,],na.rm=T)*1.1), xlab="Year", col='white', main='')
    minor.tick(nx=5, ny=2, tick.ratio=.2)
    polygon(c(F_times,rev(F_times)), c(hpd_G[2,], rev(hpd_G[1,])), col="grey90",
            border = NA)
    lines(F_times, medians_G[1,], type='l')

    # /* plot samplingProportion */
    plot(0, ylab=expression(s), xlim=c(recent-width, recent),
         ylim=c(0,max(hpd_H[2,],na.rm=T)*1.1),  xlab="Year", col='white', main='')
    minor.tick(nx=5, ny=2, tick.ratio=.2)
    polygon(c(F_times,rev(F_times)), c(hpd_H[2,], rev(hpd_H[1,])), col="grey90",
            border=NA)
    lines(F_times, medians_H[1,], type='l')


    # /* copy to pdf file, filename includes seed from logfile */
    seed <- unlist(strsplit(tail(unlist(strsplit(loglist[i,],'_')),1),'.log'))[1]
    dev.copy2pdf(file=paste("bdsky_plot_", i, "_", seed, ".pdf", sep=''))
    detach(get(paste("log", i, sep='')))
  }
}





