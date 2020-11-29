setwd("/pathToBedFile/BedFile")

# read bed file as a data frame
bed <- read.table("final_BED.txt",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

# create new data frame to store accessible chromosome coordinates, ATAC and ChIP-seq peak heights 
newBedMine <- data.frame(peaks=integer(n), startPeakPos=integer(n), endPeakPos=integer(n), a1Patac=integer(n), a1Natac=integer(n), chip1=integer(n), chip2=integer(n));

# eliminating rows from initial data frame where peak height (number of reads) of either ATAC-seq Ascl1 induced or uninduced are less than or equal to 0 and putting remaining rows 
# into newly created data frame newBedMine
newBedMine <- bed[bed$V4>0&bed$V5>0,]

# visualizing the new data frame
newBedMine

# Following steps are to normalize the total number of reads in the two columns ‘ATAC-seq peaks with Ascl1 induction’ and ‘ATAC-seq peaks without Ascl1 induction’

# Determining the sum of peak heights (sum of number of reads) for ‘ATAC-seq peaks with Ascl1 induction’, for each accessible chromatin region
a1PatacSum<- sum(newBedMine$V4)

# Visualising the sum of peak heights for ‘ATAC-seq peaks with Ascl1 induction’
a1PatacSum

# Determining the sum of peak heights (sum of number of reads) for ‘ATAC-seq peaks without Ascl1 induction’, for each accessible chromatin region
a1NatacSum<- sum(newBedMine$V5)

# Visualising the sum of peak heights for ‘ATAC-seq peaks without Ascl1 induction’
a1NatacSum

# Determining ratio of ‘total number of reads for ATAC-seq peaks without Ascl1 induction’ by ‘ total number of reads for ATAC-seq peaks with Ascl1 induction’. 
# This ratio will be used for the normalisation.
AtacPeakRatio <- a1NatacSum/a1PatacSum

# Visualising the normalisation ratio obtained.
AtacPeakRatio

# Creating new column for ‘ATAC-seq peaks without Ascl1 induction’ multiplied by ‘normalisation ratio’ since the 
# ‘total number of reads for ATAC-seq peaks without Ascl1 induction’ was greater than
# ‘total number of reads for ATAC-seq peaks with Ascl1 induction’ by a factor of ratio. And rounding the values obtained in this new column.
newBedMine$a1PProdRatio <- round(newBedMine$V4*AtacPeakRatio)

# Typecasting values in the new column as an integer
newBedMine$a1PProdRatio <- intéger(newBedMine$a1PProdRatio)

# Visualizing contents of the new column
newBedMine$a1PProdRatio

# Appending new column to the newBedMine data frame
newBedMine <- newBedMine[, c("V1", "V2", "V3", "V4", "V5", "a1PProdRatio", "V6", "V7)]

# Visualising the data frame with the appended column
newBedMine

# Creating a new data frame newBedMineFinal containing filtered data and the ratio "ATAC-seq induced" vs "uninduced peak heights" column
newBedMineFinal <- data.frame(peaks=integer(n), startPeakPos=integer(n), endPeakPos=integer(n), a1Patac=integer(n), a1Natac=integer(n), a1PProdRatio=integer(n), 
chip1=integer(n), chip2=integer(n));


# Extracting only those accessible chromatin regions (peaks) where ‘ATAC-seq peaks with Ascl1 induction’ and ‘ATAC-seq peaks without Ascl1 induction’ 
# are greater than 10 as a filter
newBedMineFinal <- newBedMine[newBedMine$V5>10&newBedMine$a1PProdRatio>10,]

# Visualising the updated data frame
newBedMineFinal

# Determing the ratio "ATAC-seq induced" vs "uninduced peak heights" (number of reads), for each peak upto 3 decimal places. This ratio is put into a new column
and is a measure of the change in chromatin accessibility, with Ascl1(induced) versus without Ascl1(uninduced).
newBedMineFinal$ratio <- signif(newBedMineFinal$a1PProdRatio / newBedMineFinal$V5,digits = 3)

# Visualising the ratio "ATAC-seq induced" vs "uninduced peak heights" column
newBedMineFinal$ratio

# Appending the ratio "ATAC-seq induced" vs "uninduced peak heights" column to the dataframe
newBedMineFinal <- newBedMineFinal[, c("V1", "V2", "V3", "V4", "V5", "a1PProdRatio", "ratio", "V6", "V7")]

# Visualising the updated data frame
newBedMineFinal

# Sorting the data frame in descending order of ratio to obtain "fold increase in chromatin accessibility, Ascl1 induced vs uninduced" 
# and putting into a new data frame
newBedMine_sorted_ratio_desc <- newBedMineFinal[with(newBedMineFinal, order(-ratio)), ]

# Visualising contents of the new data frame
newBedMine_sorted_ratio_desc 


# Following lines of code were to apply a filter so that either only those peaks are accepted where peak heights of "ATAC-seq peaks without Ascl1 induction" 
# is greater than 20 or normalized peak heights of "ATAC-seq peaks with Ascl1 induction" is greater than 20
# and put into a new data frame namely newBedMine_sorted_ratio_desc_filtered

newBedMine_sorted_ratio_desc_filtered <-data.frame(peaks=integer(n), startPeakPos=integer(n), endPeakPos=integer(n), a1Patac=integer(n), a1Natac=integer(n),
a1PProdRatio=integer(n), chip1=integer(n), chip2=integer(n));
newBedMine_sorted_ratio_desc_filtered <- newBedMine_sorted_ratio_desc[newBedMine_sorted_ratio_desc$V5>20|newBedMine_sorted_ratio_desc$a1PProdRatio>20,]
newBedMine_sorted_ratio_desc_filtered


# Following lines of code were to apply a filter to remove those peaks where peak heights are more than 100 for "ChIP-seq control" and "ChIP-seq induced" 
# and put into a new data frame namely newBedMineChipEliminated

newBedMineChipEliminated <- data.frame(peaks=integer(n), startPeakPos=integer(n), endPeakPos=integer(n), a1Patac=integer(n), a1Natac=integer(n), 
a1PProdRatio=integer(n), chip1=integer(n), chip2=integer(n));
newBedMineChipEliminated <- newBedMine_sorted_ratio_desc_filtered[newBedMine_sorted_ratio_desc_filtered$V6<100|newBedMine_sorted_ratio_desc_filtered$V7<100,]
newBedMineChipEliminated


# Generation of subset of only the first 100 peaks for barplot generation
newBedMineChipEliminated_subset<- newBedMineChipEliminated[1:100,]
newBedMineChipEliminated_subset
write.table(newBedMineChipEliminated_subset, file="readCountSubset.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE )

# Code for barplot generation
# Barplot Generation for "Fold increase in chromatin accessibility"
par(mar=c(4, 4, 1, 1))
xLabLocs1 <- barplot(newBedMineChipEliminated_subset$ratio,ylim = c(0,50), xlim=c(0, 180),width =
1.54,yaxt="none", xlab = "Peaks", ylab = "Peak Heights")
lab<-c("1","10","20","30","40","50","60","70","80","90","100")
axis(1, at=xLabLocs[c(1,10,20,30,40,50,60,70,80,90,100)],labels=lab)
axis(2,las=2)
dev.copy(jpeg,filename="BarPlotRatio.jpg");
dev.off (); 


# Barplot Generation for ATAC-seq ASCL1 uninduced
par(mar=c(4, 4, 1, 1))
xLabLocs2 <- barplot(newBedMineChipEliminated_subset$V5,ylim = c(0,900), xlim=c(0, 180),width =
1.54,yaxt="none", xlab = "Peaks", ylab = "Peak Heights")
lab<-c("1","10","20","30","40","50","60","70","80","90","100")
axis(1, at=xLabLocs[c(1,10,20,30,40,50,60,70,80,90,100)],labels=lab)
axis(2,las=2)
dev.copy(jpeg,filename="BarPlotAscl1ATAC- .jpg");
dev.off (); 


# Barplot Generation for ATAC-seq ASCL1 induced
par(mar=c(4, 4, 2, 1))
xLabLocs3 <- barplot(newBedMineChipEliminated_subset$a1PProdRatio,ylim = c(0,800), xlim=c(0, 180),width =
1.54,yaxt="none", xlab = "Peaks", ylab = "Peak Heights")
lab<-c("1","10","20","30","40","50","60","70","80","90","100")
axis(1, at=xLabLocs[c(1,10,20,30,40,50,60,70,80,90,100)],labels=lab)
axis(2,las=2)
dev.copy(jpeg,filename="BarPlotAscl1ATAC+.jpg");
dev.off (); 


# Barplot Generation for ChIP-seq 
par(mar=c(4, 4, 1, 1))
xLabLocs4 <- barplot(newBedMineChipEliminated_subset$V11, xlim=c(0, 180),width = 1.54,yaxt="none", xlab =
"Peaks", ylab = "Peak Heights")
lab<-c("1","10","20","30","40","50","60","70","80","90","100")
axis(1, at=xLabLocs[c(1,10,20,30,40,50,60,70,80,90,100)],labels=lab)
axis(2,las=2)
dev.copy(png,filename="BarPlotChIP.png");
dev.off (); 


# Generation of subset of only the first 1000 peaks for usage in enrichment analysis
newBedMineChipEliminated_subset_enrichr<- newBedMineChipEliminated[1:1000,]
newBedMineChipEliminated_subset_enrichr


# Command to write BED file having 1000 peaks for enrichment analysis into a txt file
write.table(newBedMineChipEliminated_subset_enrichr, file="BEDEnrichr.txt", quote=FALSE, sep="\t", row.names =
FALSE, col.names = FALSE )


# Command to write final BED file having all peaks sorted in descending order of ratio into a txt file
write.table(newBedMineChipEliminated, file="BEDFoldIncrease.txt", quote=FALSE, sep="\t", row.names = FALSE,
col.names = FALSE )





















