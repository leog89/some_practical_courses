# Exercice 2 Feuille TD2

# Calcul facteur de Bayes

# On reprend l'exemple du modèle linéaire multivarié
# données dans R "LifeCycleSavings" 
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

#
matinvXX=solve(t(X)%*%X)


# dans l'expression de p(Y|X) il suffit de conserver le dénominateur (g+1)^(p/2) * [...]^{a+n/2} 
# (tout le reste s'annule dans le Facteur de Bayes)

# pour modèle total  (attention on calcule log de l'expression)
pourTout=log( (b+0.5*(sum(y*y)- (g/(g+1))*(t(y)%*%X%*%matinvXX%*%t(X)%*%y)[1,1]  ) )^(a+0.5*n) )
pourTout=pourTout+(0.5*dim(X)[2])*log(g+1)

# pour modèle sans pop15
Xpop15=X[,-2]
matinvXXpop15=solve(t(Xpop15)%*%Xpop15)
pourSansPop15=log( (b+0.5*(sum(y*y)- (g/(g+1))*(t(y)%*%Xpop15%*%matinvXXpop15%*%t(Xpop15)%*%y)[1,1]  ) )^(a+0.5*n)  )
pourSansPop15=pourSansPop15+(0.5*dim(Xpop15)[2])*log(g+1)

# pour modèle sans pop75
Xpop75=X[,-3]
matinvXXpop75=solve(t(Xpop75)%*%Xpop75)
pourSansPop75=log( (b+0.5*(sum(y*y)- (g/(g+1))*(t(y)%*%Xpop75%*%matinvXXpop75%*%t(Xpop75)%*%y)[1,1]  ) )^(a+0.5*n)  )
pourSansPop75=pourSansPop75+(0.5*dim(Xpop75)[2])*log(g+1)

# pour modèle sans dpi
Xdpi=X[,-4]
matinvXXdpi=solve(t(Xdpi)%*%Xdpi)
pourSansdpi=log( (b+0.5*(sum(y*y)- (g/(g+1))*(t(y)%*%Xdpi%*%matinvXXdpi%*%t(Xdpi)%*%y)[1,1]  ) )^(a+0.5*n)  )
pourSansdpi=pourSansdpi+(0.5*dim(Xdpi)[2])*log(g+1)

# pour modèle sans ddpi
Xddpi=X[,-5]
matinvXXddpi=solve(t(Xddpi)%*%Xddpi)
pourSansddpi=log( (b+0.5*(sum(y*y)- (g/(g+1))*(t(y)%*%Xddpi%*%matinvXXddpi%*%t(Xddpi)%*%y)[1,1]  ) )^(a+0.5*n)  )
pourSansddpi=pourSansddpi+(0.5*dim(Xddpi)[2])*log(g+1)

# Pour  rappel résultats fonction lm :
print(summary(reslm))

# Présentation log BF
cat("log BF pour test de Beta=0","\n")
cat(" pop15 :",pourTout-pourSansPop15,"\n")
cat(" pop75 :",pourTout-pourSansPop75,"\n")
cat(" dpi :",pourTout-pourSansdpi,"\n")
cat(" ddpi :",pourTout-pourSansddpi,"\n")


##################################################################################################
# Alternative : utilisation library BAS (Bayesian Variable Selection and Model Averaging using Bayesian Adaptive Sampling)
library(BAS)

resBAS=bas.lm(formula = sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings,prior="g-prior",alpha=50)
print(resBAS)
print(summary(resBAS))

# si beaucoup de covariables utilise MCMC
