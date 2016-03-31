#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 16/05/2013 by Joe Gallagher 
# joedgallagher@gmail.com
#
# This script contains a function to calculate the maximal rate of
# the linear phase of a reaction, Vmax, from optical density data, 
# and a diagnostic function to visualise smoothed data and the 
# fitted Vmax gradient.
#
# Examples are provided at the end of the document.
#
# This script requires a dataframe with defined "OD" and 
# "time" columns, such as that derived from a raw Softmax V-5.1 
# ouput file, which can itself be processed with the use of another script, # "SoftmaxToRX.R", made by Quentin Geissmann (qgeissmann@gmail.com).
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

require(chron)
graphics.off()

source("SoftmaxToR9.R")

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
calcVmaxIntervals <- function(data, windowSize = 10, smoothing = 31, timeStart = 0, timeEnd = globalXmax){
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Function performs a running median to smooth data, before
# bootstrapping the dataset into overlapping windows of
# length windowSize, and finding the window with the largest linear
# rate of reaction, Vmax.
# Outputs a dataframe containing Vmax, row and column.
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# data = dataframe outputted by SofToR() function
# windowSize = window to calculate linear rate of reaction over
# smoothing = size of running median window to smooth raw data
# timeStart, timeEnd = times to subset data by
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
	allVmaxs <- matrix(numeric(0), ncol = 3)
	nrows <- length(unique(data$row))
	ncols <- length(unique(data$column))
	globalXmax <- max(data$time)
	
	#loop through wells
	for(i in 1:nrows) for(j in 1:ncols){
		#subset
		dd <- subset(data, row==i & column==j & time>=timeStart & time<=timeEnd)	
		#time
		x <- dd$time
		time.int <- median(diff(dd$time))
		#smooth OD
		y <- runmed(dd$OD, smoothing)

		#create lists of individual windows
		df <- data.frame(x,y)
		start <- seq(1, nrow(df), by = 1)
		list <- lapply(start, function(i) df[c(i:(i+windowSize-1)),])
		#omit windows containing NAs
		list <- list[sapply(list, function(i) nrow(na.omit(i)) == windowSize)]
		#find window with maximum rate of OD change
		VmaxElement <- list[[which.max(lapply(list, function(z) z$y[windowSize] - z$y[1]))]]
		#calculate Vmax (change in OD per minute)
		Vmax <- (VmaxElement$y[nrow(VmaxElement)] - VmaxElement$y[1]) / (VmaxElement$x[nrow(VmaxElement)] - VmaxElement$x[1]) * 60
		
		#append
		dfi <- data.frame(Vmax, row = i, col = j)
		allVmaxs <- rbind(allVmaxs, dfi)
	}
	return(allVmaxs)
}

calcVmaxRegression <- function(data, smoothing = 31, timeStart = 0, timeEnd = globalXmax, filename = "Rplots.pdf"){
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Function simply performs a linear regression on OD data to determine 
# rate of linear change, and outputs a dataframe containing Vmax, row 
# and column.
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# data = dataframe outputted by SofToR() function
# windowSize = window to calculate Vmax over
# smoothing = running median window to apply to raw data
# timeStart, timeEnd = times to subset data by
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
	allVmaxs <- matrix(numeric(0), ncol = 3)
	nrows <- length(unique(data$row))
	ncols <- length(unique(data$column))
	globalXmax <- max(data$time)
	
	#loop through wells
	for(i in 1:nrows) for(j in 1:ncols){
		#subset
		dd <- subset(data, row==i & column==j & time>=timeStart & time<=timeEnd)
		x <- dd$time
		y <- runmed(dd$OD, smoothing)
		time.int <- median(diff(dd$time))
		
		#calculate Vmax
		modVmax <- lm(y ~ x)
		Vmax <- modVmax$coefficients[2] * 60
		
		#append
		dfi <- data.frame(Vmax, row = i, col = j)
		allVmaxs <- rbind(allVmaxs, dfi)
	}
	return(allVmaxs)
}

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
plotVmaxIntervals <- function(data, windowSize = 10, smoothing = 31, timeStart = 0, timeEnd = globalXmax, filename = "Rplots.pdf"){
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Function plots smoothed OD data and fitted Vmax (calculated by 
# maximal interval linear change) for each well, and outputs a pdf
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# data = dataframe outputted by SofToR() function
# windowSize = window to calculate Vmax over
# smoothing = running median window to apply to raw data
# timeStart, timeEnd = times to subset data by
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
	pdf(filename)
	nrows <- length(unique(data$row))
	ncols <- length(unique(data$column))
	globalYmax <- max(data$OD)
	globalXmax <- max(data$time)
	
	#loop through wells
	for(i in 1:nrows) for(j in 1:ncols){
		#subset
		dd <- subset(data, row==i & column==j & time>=timeStart & time<=timeEnd)
		#time
		x <- dd$time
		time.int <- median(diff(dd$time))
		#smooth OD
		y <- runmed(dd$OD, smoothing)
		
		df <- data.frame(x,y)
		start <- seq(1, nrow(df), by = 1)
		
		#make list of individual windows
		list <- lapply(start, function(i) df[c(i:(i+windowSize-1)),])
		#omit lists containing NAs
		list <- list[sapply(list, function(i) nrow(na.omit(i)) == windowSize)]
		
		#find window with maximum rate of OD change
		VmaxElement <- list[[which.max(lapply(list, function(z) z$y[windowSize] - z$y[1] ))]]
		#calculate Vmax (change in OD per minute)
		Vmax <- (VmaxElement$y[nrow(VmaxElement)] - VmaxElement$y[1]) / (VmaxElement$x[nrow(VmaxElement)] - VmaxElement$x[1]) * 60
		
		#make plots
		plot(y ~ x,
			cex=0.5, pch=19,
			xlab = "Time (s)", ylab = "OD",
			main=sprintf("row = %i\ncol = %i",i,j),
			ylim=c(0,globalYmax), xlim=c(0,globalXmax),
			type= 'o'
		)
		abline(lm(VmaxElement$y ~ VmaxElement$x), lwd = 2, col = "red")
		#report Vmax per minute
		text(globalXmax, 0, paste("Vmax = ", signif(Vmax, 4)), pos=2)
	}
	dev.off()
}

plotVmaxRegression <- function(data, smoothing = 31, timeStart = 0, timeEnd = globalXmax, filename = "Rplots.pdf"){
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Function plots smoothed OD data and fitted Vmax (calculated by
# linear regression) for each well, and outputs a pdf
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# data = dataframe outputted by SofToR() function
# windowSize = window to calculate Vmax over
# smoothing = running median window to apply to raw data
# timeStart, timeEnd = times to subset data by
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
	pdf(filename)
	nrows <- length(unique(data$row))
	ncols <- length(unique(data$column))
	globalYmax <- max(data$OD)
	globalXmax <- max(data$time)
	
	#loop through wells
	for(i in 1:nrows) for(j in 1:ncols){
		#subset
		dd <- subset(data, row==i & column==j & time>=timeStart & time<=timeEnd)
		x <- dd$time
		y <- runmed(dd$OD, smoothing)
		time.int <- median(diff(dd$time))
		
		#calculate Vmax (change in OD per second, for plotting)
		modVmax <- lm(y ~ x)
		Vmax <- modVmax$coefficients[2]
		
 		#make plots
		plot(y ~ x, cex=0.5, pch=19, xlab = "Time (s)", ylab = "OD", main=sprintf("row = %i\ncol = %i",i,j), ylim=c(0,globalYmax), xlim=c(0,globalXmax), type= 'o')
		abline(modVmax, lwd = 2, col = "red")
		#report Vmax per minute
		text(globalXmax, 0, paste("Vmax = ", signif(Vmax*60, 4)), pos=2)
	}
	dev.off()
}

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
## EXAMPLES
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# #Declare functions
#source("VmaxFunctions.R")

# #Read data
# d <- sofToR("example.text")

# #Pre-subset data by wells, e.g.
# d = subset(d, d$row==1 & d$column%in%1:4)

# #Pre-subset data by time, e.g.
# d = subset(d, d$time >= 2000)

# #Calculate a linear regression Vmax (e.g. window size of 10 reads, smoothing window size of 31 reads, data subsetted by time)
#calcVmaxIntervals(d, windowSize = 10, smoothing = 31, timeStart = 1000, timeEnd = 2000)

# #...and visualise
#plotVmaxIntervals(d, windowSize = 10, smoothing = 31, timeStart = 1000, timeEnd = 2000, filename = "VmaxIntervalExample.pdf")

# #Calculate a linear regression Vmax (e.g. window size of 20 reads, data subsetted by time)
#calcVmaxRegression(d, windowSize = 20, timeStart = 2000)

# #...and visualise
#plotVmaxRegression(d, windowSize = 20, timeStart = 2000, filename = "VmaxRegressionExample.pdf")

# #Write Vmax dataframe to csv
# write.csv(calcVmax(d), file="exampleVmax.csv")

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
########################### END ##############################