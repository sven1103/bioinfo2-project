
library(MASS)
library(ggplot2)
library(gplots)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(9,'YlOrRd')))
r <- rf(20)
torsion_helices <- read.table("torsion_angles_helices.csv", sep=",")
torsion_sheets <- read.table("torsion_angles_sheets.csv", sep=",")
torsion_all <- read.table("torsion_angles_all.csv", sep=",")

h2_helices <- hist2d(torsion_helices[,1], torsion_helices[,2], nbins=100, col=r)
h2_sheets <- hist2d(torsion_sheets[,1], torsion_sheets[,2], nbins=100, col=r)
h2_all <- hist2d(torsion_all[,1], torsion_all[,2], nbins=100, col=r)
density_helices <- kde2d(torsion_helices[,1], torsion_helices[,2], h=20, n=100, lims = c(range(torsion_helices[,1]), range(torsion_helices[,2])))
density_sheets <- kde2d(torsion_sheets[,1], torsion_sheets[,2], h=20, n=100, lims = c(range(torsion_sheets[,1]), range(torsion_sheets[,2])))
density_all <- kde2d(torsion_all[,1], torsion_all[,2], h=20, n=100, lims = c(range(torsion_all[,1]), range(torsion_all[,2])))

filled.contour(h2_helices, z=density_helices$z, xlab = "phi torsion angles", ylab = "psi torsion angles",
                 levels = pretty(c(0.0000001, 0.00001), 20), col=r, main="Oberserved torsion values in alpha helices")
filled.contour(h2_sheets, z=density_sheets$z, xlab = "phi torsion angles", ylab = "psi torsion angles",
               levels = pretty(c(0.0000001, 0.00001), 20), col=r, main="Oberserved torsion values in sheets")
filled.contour(h2_all, z=density_all$z, xlab = "phi torsion angles", ylab = "psi torsion angles",
               levels = pretty(c(0.0000001, 0.00001), 20), col=r, main="Oberserved torsion values in all residues")

plot(torsion_all[,1], torsion_all[,2], pch=21, main="Oberserved torsion values in all residues")
plot(torsion_sheets[,1], torsion_sheets[,2], pch=21, main="Oberserved torsion values in all residues")
plot(torsion_helices[,1], torsion_helices[,2], pch=21, main="Oberserved torsion values in all residues")
