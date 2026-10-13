
# simulation échantillon loi de Cauchy mu=2
n =5 # taille échantillon
x=rcauchy(n,location=2)

# grille de valeurs pour mu 
valmu=seq(from=(min(x)-2),to=(max(x)+2),length.out=100)

# calcul vraisemblance pour grille de valeurs
vraisemblance<-rep(1,length(valmu))
for (i in 1:n) vraisemblance=vraisemblance*dcauchy(valmu,location=x[i])

# fonctionvraisemblance
foncvrais<-function(valeur){
	vrais=1
	for (i in 1:n) vrais=vrais*dcauchy(valeur,location=x[i])
	vrais
}

# on standardise (travail à une constante multiplicative près, pour plot)
vraisemblance=vraisemblance/max(vraisemblance)

# priori Gaussienne N(0,10)
densitépriori=dnorm(valmu,mean=0,sd=10^0.5)
densitépriori=densitépriori/max(densitépriori)

# Loi a posteriori (toujours à une constante multiplicative près)
densitéposteriori= vraisemblance*densitépriori
densitéposteriori= densitéposteriori/max(densitéposteriori)

# plot
#plot(valmu,vraisemblance,type="l")
plot(valmu,densitéposteriori,ylab=" ",type="l",yaxt="n",lwd=2)
lines(valmu,0.3*densitépriori,lty=2)
points(x,rep(0,n),pch=19)


# Estimation moyenne loi a posteriori par Importance Sampling

# On utilise comme loi intrumentale la loi a priori (choix naif)
nbresimul=100
simul=rnorm(nbresimul,mean=0,sd=10^0.5)
# calul des poids
w=rep(0,nbresimul)
for (k in 1:nbresimul) w[k]=foncvrais(simul[k])
w=w/dnorm(simul,mean=0,sd=10^0.5)
# on standardise les poids
w=w/sum(w)
points(simul,rep(0,nbresimul),col="red")
points(simul,w/max(w),col="red",type="h")

moy=sum(w*simul)
cat("Estimation moyenne loi a posteriori :",moy)

lines(c(moy,moy),c(0,1000),col="blue",lwd=3)
legend("topleft",legend=c("priori/instru.","posteriori","simulations","estimation moyenne"),col=c("black","black","red","blue"),lty=c(2,1,1,1),bty="n")


