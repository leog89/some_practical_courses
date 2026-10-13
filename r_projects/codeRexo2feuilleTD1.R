
# Exercice 2 Feuille TD1 Application numérique

######################
######################
#  Modele linéaire simple  
######################
######################


# simulation observations Y=beta_0+beta_1*x+epsilon
beta_0=0
beta_1=1
n=30
x=runif(n,0,1)
sigma=0.9
set.seed(2025)
y=beta_0+beta_1*x+rnorm(n,mean=0,sd=sigma)

X=cbind(rep(1,n),x)


###############################
# I - On suppose d'abord sigma connu

# résultats fonction lm
reslm=lm(y~x)
estibeta_0=reslm$coefficients[1]
estibeta_1=reslm$coefficients[2]

# Bayésien avec g-prior (g à fixer, g=n par défaut)
g=n  
# loi a posteriori pour beta : N_2(mu,Sigma)
# avec 
mu=(g/(g+1))*solve(t(X)%*%X)%*%t(X)%*%y
Sigma=(sigma^2)*(g/(g+1))*solve(t(X)%*%X)


# plot résultat estimation
plot(x,y,pch=19)
sortx=sort(x)
lines(sortx,beta_0+beta_1*sortx,lty=2)
lines(sortx,estibeta_0+estibeta_1*sortx,col="blue",lwd=2)
lines(sortx,mu[1]+mu[2]*sortx,col="red",lwd=2)


# intervalles de confiance fréq (à 95%)
# à partir de la loi de l'estimateur moindres carrés
# pour beta_0 loi N(estibeta_0,sigma^2*(X'X)^{-1}_(1,1))
matinvXX=solve(t(X)%*%X)
IC0=c(qnorm(0.025,mean=estibeta_0,sd=sigma*(matinvXX[1,1])^0.5), qnorm(0.975,mean=estibeta_0,sd=sigma*(matinvXX[1,1])^0.5))
IC1=c(qnorm(0.025,mean=estibeta_1,sd=sigma*(matinvXX[2,2])^0.5), qnorm(0.975,mean=estibeta_1,sd=sigma*(matinvXX[2,2])^0.5))

# intervalles de crédibilité bayésien (à 95%)
CI0=c(qnorm(0.025,mean=mu[1],sd=(Sigma[1,1])^0.5), qnorm(0.975,mean=mu[1],sd=(Sigma[1,1])^0.5))
CI1=c(qnorm(0.025,mean=mu[2],sd=(Sigma[2,2])^0.5), qnorm(0.975,mean=mu[2],sd=(Sigma[2,2])^0.5))

cat("\n", "\n \n \n", "********************************************* \n ********************************************* \n Exemple de régression linéaire simple  \n", "données  simulées (p=2)  \n *********************************************\n *********************************************\n", "\n")
cat(" * Cas sigma^2 connu \n","\n")
# présentation résultats estimation
cat(" Estimateur mondres-carrés \n"," de beta_0 : ",estibeta_0,"\n  de beta_1 : ",estibeta_1,"\n Estimateur bayésien    ( g =",g,")\n"," de beta_0 : ",mu[1],"\n  de beta_1 : ",mu[2])
# présentation résultats intervalles confiance/crédibilité
cat("\n","\n")
cat(" Intervalle de confiance à 95%\n"," pour beta_0: ",IC0[1],IC0[2],"\n  pour beta_1: ",IC1[1],IC1[2],"\n Intervalle de crédibilité à 95%   ( g =",g,")\n"," pour beta_0: ",CI0[1],CI0[2],"\n  pour beta_1: ",CI1[1],CI1[2])


# région de confiance
# utilisation library "ellipse"
library(ellipse)
# la fonction ellipse donne ellipsoide de confiance pour une gaussienne (par défaut à 95%)
resMLE=ellipse(sigma*matinvXX,center=c(estibeta_0,estibeta_1))
resBayesien=ellipse(Sigma,center=mu)
resPrior=ellipse(g*sigma^2*matinvXX,center=c(0,0))

# plot des résultats

# comparaison MLE/moindres-carrés avec Bayésien
minx=min(c(resMLE[,1],resBayesien[,1]))
maxx=max(c(resMLE[,1],resBayesien[,1]))
miny=min(c(resMLE[,2],resBayesien[,2]))
maxy=max(c(resMLE[,2],resBayesien[,2]))
# plot
quartz()                               # Nouvelle fenetre graphique (pour MAC)
plot(resMLE[,1],resMLE[,2],type="l",col="blue",lwd=2,xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="beta_0",ylab="beta_1")
lines(resBayesien[,1],resBayesien[,2],col="red",lwd=2)
lines(c(-100,100),c(0,0))
lines(c(0,0),c(-100,100))
points(beta_0,beta_1,pch=4)
points(estibeta_0,estibeta_1,col="blue",pch=4)
points(mu[1],mu[2],col="red",pch=4)
legend("bottomleft",legend=c("MLE","Bayes"),col=c("blue","red"),lty=c(1,1))


# comparaison a priori / a posteriori
minx=min(c(resPrior[,1],resBayesien[,1]))
maxx=max(c(resPrior[,1],resBayesien[,1]))
miny=min(c(resPrior[,2],resBayesien[,2]))
maxy=max(c(resPrior[,2],resBayesien[,2]))
# plot
quartz()                               # Nouvelle fenetre graphique (pour MAC)
plot(resPrior[,1],resPrior[,2],type="l",col="red",lty=2,lwd=2,xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="beta_0",ylab="beta_1")
lines(resBayesien[,1],resBayesien[,2],col="red",lwd=2)
lines(c(-100,100),c(0,0))
lines(c(0,0),c(-100,100))
points(beta_0,beta_1,pch=4)
points(estibeta_0,estibeta_1,col="blue",pch=4)
points(mu[1],mu[2],col="red",pch=4)
legend("bottomleft",legend=c("prior","posterior"),col=c("red","red"),lty=2:1)

###################################
# I - On suppose maintenant sigma inconnu
# loi a priori sur sigma^2 : inverse gamma (a,b)
a=1
b=1

# loi a posteriori sur beta : loi student multivariée de vecteur mu et de matrice S et de ddl nu
nu=2*a+n
S=(2*b+sum(y*y)-(g/(g+1))*(t(y)%*%X%*%matinvXX%*%t(X)%*%y)[1,1] )*(g/(g+1))*matinvXX/nu

bornesT=c(qt(0.025,df=nu),qt(0.975,df=nu))
CIvar0=bornesT*(S[1,1])^0.5+mu[1]
CIvar1=bornesT*(S[2,2])^0.5+mu[2]

cat("\n \n","* Cas sigma^2 inconnu \n","\n")
# intervalles de confiances pour lm
cat("Intervalles de confiance à 95% \n")
print(confint(reslm))
# pour Bayesien
# présentation résultats intervalles confiance/crédibilité
cat("\n")
cat(" Intervalle de crédibilité à 95%   ( g =",g,")\n"," pour beta_0: ",CIvar0[1],CIvar0[2],"\n  pour beta_1: ",CIvar1[1],CIvar1[2])


# Question sup : faut-il éliminer beta_0 ? beta_1 ? Les deux ?

# Question sup : prédiction de y si x=0.5. Densité de prédiction ?
# On trouve loi de Student. 


##############################
##############################
#  Modèle régression linéaire multiple  
##############################
##############################

# Exemple : données dans R "LifeCycleSavings" 
# (exemple utilisé dans des aides pour fonction lm)

n=50

# variable réponse : variable sr
y=LifeCycleSavings$sr

# Les 4 covariables : pop15 + pop75 + dpi + ddpi
X=cbind(rep(1,n),LifeCycleSavings$pop15,LifeCycleSavings$pop75,LifeCycleSavings$dpi,LifeCycleSavings$ddpi)

# ici p=5

# résultat fonction lm
reslm=lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
estibeta=reslm$coefficients

# Bayésien avec g-prior (g à fixer, g=n par défaut)
g=n  
# loi a priori sur sigma^2 : inverse gamma (a,b)
a=1
b=1

# loi a posteriori pour beta : Student centrée en
mu=(g/(g+1))*solve(t(X)%*%X)%*%t(X)%*%y
#  et de matrice S et de ddl nu
nu=2*a+n
matinvXX=solve(t(X)%*%X)
S=(2*b+sum(y*y)-(g/(g+1))*(t(y)%*%X%*%matinvXX%*%t(X)%*%y)[1,1] )*(g/(g+1))*matinvXX/nu

bornesT=c(qt(0.025,df=nu),qt(0.975,df=nu))
CIconst=bornesT*(S[1,1])^0.5+mu[1]
CIpop15=bornesT*(S[2,2])^0.5+mu[2]
CIpop75=bornesT*(S[3,3])^0.5+mu[3]
CIdpi=bornesT*(S[4,4])^0.5+mu[4]
CIddpi=bornesT*(S[5,5])^0.5+mu[5]

cat("\n", "\n \n \n", "********************************************* \n ********************************************* \n Exemple de régression linéaire multivariée  \n", "données  LifeCycleSavings  \n *********************************************\n *********************************************\n", "\n")

cat("Estimation \n")
# présention Estimation
print(summary(reslm))
cat("Estimateur bayésien :",mu)

cat("\n","\n")
# intervalles de confiances pour lm
cat("intervalles de confiance à 95%\n")
print(confint(reslm))
# pour Bayesien
# présentation résultats intervalles confiance/crédibilité
cat("\n","Intervalle de crédibilité à 95%   ( g =",g,")\n"," pour const: ",CIconst[1],CIconst[2],"\n  pour pop15: ",CIpop15[1],CIpop15[2],"\n  pour pop75: ",CIpop75[1],CIpop75[2],"\n  pour dpi: ",CIdpi[1],CIdpi[2],"\n  pour ddpi: ",CIddpi[1],CIddpi[2])

