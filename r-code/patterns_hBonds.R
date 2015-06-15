
hbonds <- read.table("1q48.hb", sep=" ")

filteredBonds <- hbonds[hbonds$V3 < -1 | hbonds$V4 < -1,]

filteredBondsColoring <- unlist(apply(filteredBonds, 1, colorAssigner))


plot(filteredBonds$V1 ~ filteredBonds$V2, main="Hydrogen Bonds Pattern encode Secondary Structures", xlab="primary Sequence",
     ylab="primary Sequence")
abline(a=0, b=1, col="blue")
abline(a=-4, b=1, col="red")