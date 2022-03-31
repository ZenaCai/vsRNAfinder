suppressMessages(library(pracma))
options(scipen = 200)
args = commandArgs(T)
filename <- args[1]
name <- args[2]
strand <- args[3]
peakheight <- args[4]
setwd(name)

data <- read.table(filename, sep = '\t', header = TRUE)
pp = "[+]{1,}[0]*[-]{1,}"

peak <- findpeaks(data$y_smooth, peakpat=pp, minpeakheight = as.numeric(peakheight))
data_peak <- cbind(peak, data[peak[, 2],], data[peak[, 3],]$X1, data[peak[, 4],]$X1)
outfile <- paste('peak.y_smooth.', strand, '.txt', sep = '')
write.table(data_peak, outfile, 
            sep = '\t', row.names = FALSE, col.names = FALSE)
