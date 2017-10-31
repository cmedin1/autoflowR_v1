


rm(list=ls())
library(flowCore)
library(flowDensity)
library(GEOmap)
library(flowViz)
setwd("/Users/i-cubed/Dropbox (iCubed)/WRAIR/WRAIR_flow")

#looking at one file
?read.FCS
f<-read.FCS('./2136_Th17-Treg_Panel_032817_HF/data_file/Subject 2136-002_Visit 1_C03.fcs')
f
#number of events
nrow(f)
#channel names
colnames(f)
#extract expression values into a matrix, the numbers represent measured intensities of each event/cell
E<-exprs(f)
#the dimensions of the data
dim(E)
#explore the meta data in the fcs file
f@description
names(f@description)
f@description$'TUBE NAME'
f@parameters@data
f@parameters@data[1, c("minRange","maxRange")]
#simple plot
plot(f,c("FSC-A","SSC-A"))
#without smoothing
plot(f, c("FSC-A","SSC-A"), ylim=c(0,5000), smooth=FALSE)

#to determine if each data is to be viewed on a linear 'LIN' or 'LOG', colname(f)[3] is SSC-A data
colnames(f)[3]
f@description$'P3DISPLAY'# this indicates LIN=linear



fs<-read.flowSet(file=NULL, "./2136_Th17-Treg_Panel_032817_HF/data_file", pattern=".fcs")
fs
# You can see sample names as well as the channel names
sampleNames(fs)
length(fs)
colnames(fs)
# A flowSet object is similar to a list, a list of flowFrames
fs[["Subject 2136-002_Visit 1_C03.fcs"]]
fs[[1]]
# Use fsApply to get cell counts for all samples
nrow(fs[[1]])
fsApply(fs, nrow)
# Use fsApply to extract the TUBE NAME keyword in all samples
fsApply(fs, function(f) f@description$`TUBE NAME`)
### Plotting excercise ##############
plot(fs[[2]], c("FSC-A", "SSC-A"), ylim = c(0, 5000), smooth = FALSE)
# Plot the density of the forward scatter area values for the first sample:
E <- exprs(fs[[1]])
fscValues <- E[, "FSC-A"]
fscValues[1:10]
plot(density(fscValues))
# We can plot all 3 samples on one plot:
par (mfrow = c(3, 1)) # This creates a plot region with a single column of 3 subplots
plot(fs[[1]], c("FSC-A", "SSC-A"), main = sampleNames(fs)[1], ylim = c(0, 5000), smooth=FALSE)
plot(fs[[2]], c("FSC-A", "SSC-A"), main = sampleNames(fs)[2], ylim = c(0, 5000), smooth=FALSE)
plot(fs[[3]], c("FSC-A", "SSC-A"), main = sampleNames(fs)[3], ylim = c(0, 5000), smooth=FALSE)

#compensation
?compensate

C03<- read.FCS('./2136_Th17-Treg_Panel_032817_HF/data_file/Subject 2136-002_Visit 1_C03.fcs')
C03
C03c <- C03@description$`SPILL`
C03c
C03comp <- compensate(C03, C03c)
summary(C03c)
summary(C03comp)

C04<-read.FCS('./2136_Th17-Treg_Panel_032817_HF/data_file/Subject 2136-002_Visit 2_C04.fcs')
C04
C04c <- C03@description$`SPILL`
C04c
C04comp <- compensate(C04, C04c)
summary(C04c)
summary(C04comp)






#gating

cd3 <- "APC-Cy7 R780"
dump <- "BV510 V525"
source("/Users/i-cubed/Dropbox (iCubed)/WRAIR/WRAIR_flow/support_functions.R")
pooled.frame <- getGlobalFrame(fs)
{plotDens(pooled.frame, c(cd3, dump))
  abline(v = 1, lwd=2, col = "blue")
  abline(h = 1.5, lwd=2, col="blue")}




