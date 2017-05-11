library(RADami)

> loc.dat <- read.pyRAD('/Users/fgrewe/odrive/ACD/projects/@Steve-lichen_rad/matrix_plot/c90d6m4p3_ex5000red_nobracket_nioinfo.loci')
print.pyRAD.loci(loc.dat)
# [1] "pyRAD.loci object read from /Users/fgrewe/odrive/ACD/projects/@Steve-lichen_rad/matrix_plot/c90d6m4p3_ex5000red_nobracket_nioinfo.loci on Sun Apr 23 19:34:41 2017"
# [1] "Contains:"
# [1] "-- 47 individuals"
# [1] "-- 76371 loci"
plot.locus.dist(locus.dist(loc.dat, proportional = TRUE, upper = TRUE), ladderize(tree), trW = 1, labelsW = 1, scalar = 4)
#plot was further modified in AI