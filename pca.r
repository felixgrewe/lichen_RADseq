#PCA 50% missing data
library(ape);
library(adegenet);

#read stru file
data.str <- read.table('complex_p80_nohaydenii.str', header = FALSE);
data_temp <- data.str[,3:4997];
row.names(data_temp) <- data.str[,1];
data <- df2genind(data_temp, ploidy = 1, NA.char = '-9', pop = data.str[,2]);
X <- tab(data, NA.method = 'mean');

dapc1 <- dapc(data, data.str[,2]);
#1/3 the number of inds
mycol <- c('darkblue','darkgreen','gold');
par(mfrow=c(2,1));
scatter(dapc1, col = mycol, scree.pca = TRUE, posi.pca = 'topleft', posi.da = 'bottomleft');
compoplot(dapc1, posi='bottomright', ncol = 1, col = mycol, cleg = 0.8, cex.names = 0.4);


#do pca
pca.data <- dudi.pca(X, center = TRUE, scale = TRUE);
#plot with 50 eigen values
s.class(pca.data$li, fac=pop(data), clabel = .8, col=transp(funky(10), .8), axesel=FALSE,cstar=0,cpoint=5);
add.scatter.eig(pca.data$eig[1:10],3,1,2,ratio=.2, posi = 'topleft');

