par(family="serif",bg="white")
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
##
#the prior
plot(NA,
   xlab="Parameter Value",
   ylab="Density",
   main="Parameter Prior",
   frame.plot=F,
   xaxt="n",
   yaxt="n",
   xlim=c(-4,4),
   ylim=c(0,0.4))
x1 <- seq(-4,4,0.01)
y1 <- dnorm(x=x1,mean=0,sd=1)
polygon(y=c(y1,rep(0,length(y1))),
   x=c(x1,rev(x1)),
   col="grey",
   border=F)
axis(1,col="grey")
axis(2,col="grey")

#the conflict data
conflicts_monument <- MayaConflict_1y$Conflict / MayaConflict_1y$Monuments_adj
conflicts_monument_finite <- conflicts_monument[which(is.finite(conflicts_monument))]
hist(conflicts_monument_finite,
   border="white",
   col="darkgrey",
   main="Conflict",
   xlab="Conflicts per monument",
   xaxt="n",
   yaxt="n")
axis(1,col="grey")
axis(2,col="grey")

#the climate covariate
d18O <- MayaConflict_1y[index,7]-mean(MayaConflict_1y[index,7])
hist(d18O,
   border="white",
   col="darkgrey",
   main="d18O",
   xlab="ppm deviations",
   xaxt="n",
   yaxt="n")
axis(1,col="grey")
axis(2,col="grey")
