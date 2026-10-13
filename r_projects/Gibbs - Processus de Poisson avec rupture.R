################################################################################
###### TD 3 STAT BAYESIENNE : POINT DE RUPTURE DANS UN MODELE DE POISSON #######
################################################################################

## Modules

library(GeDS)

#

### Importation des données

donnees = as.matrix(coalMining)

annees = donnees[,1]
morts = donnees[,2]
n = 112

# On définit les vecteurs des résultats

nbsimul = 5000

# lambda : le paramètre de la P(.) avant le point de rupture
# theta : le paramètre de la P(.) après le point de rupture
# T : instant de la rupture (théoriquement entre 2 et n-1 mais entre 1 et 112 + simple)

vecteurlambda = rep(NA, nbsimul)
vecteurtheta = rep(NA, nbsimul)
vecteurT = rep(NA, nbsimul)

# valeurs des lois à priori

alamb = 0.5 ; blamb = 0.1
atheta = 0.5 ; btheta = 0.1

# valeurs de départ

vecteurlambda[1] = 10
vecteurtheta[1]=10
vecteurT[1] = 20

# On lance l'algorithme de Gibbs

for(ite in 2:nbsimul){
  
  # Simulation d'un nouveau lambda
  vecteurlambda[ite] = rgamma(1, shape = (alamb + sum(morts[1:(vecteurT[(ite-1)])])), 
                                          rate = (blamb + vecteurT[(ite-1)]))[1]
  # Simulation d'un nouveau theta
  vecteurtheta[ite] = rgamma(1, shape = (atheta + sum(morts[(vecteurT[(ite-1)]+1):n] )), 
                              rate = (btheta + n - vecteurT[(ite-1)]))[1]
  
  # Simulation d'un nouveau T
  vecteurpoids = rep(NA, n)
  for(t in 1:n){
    vect=c(rep(vecteurlambda[ite], t), rep(vecteurtheta[ite], (n-t)))
    
    uiz = (exp((-1)*vect))*(vect^morts)
    
    vecteurpoids[t] = prod(uiz)
  }
  vecteurpoids=vecteurpoids/sum(vecteurpoids)
  
  vecteurT[ite]=sample(c(1:n), size=1, prob = vecteurpoids)[1]
  
  if(ite%%200==0) cat("it",ite,"\n")
}

hist(vecteurT, breaks=25)


