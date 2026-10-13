## Case Study 3 MAT R


simul_elec <- function(n, case=c("IC", "IAC_star"), B){
  
  # case should be either IC or IAC
  case<-match.arg(case)
  
  # n should be an odd number
  if(n %% 2 != 1){stop("n should be odd!")}
  
  prop_happy_voters <- matrix(0, nrow=B, ncol=2, dimnames=list(NULL, c(paste("M1 (", case, ")"), paste("M2 (", case, ")"))))
    
  total_voters <- 3*n
  for(b in 1:B){

    p<-if(case=="IC") 1/2 else runif(3,0,1)
    
    voters<-rbinom(3, n, p)
    
    py_voters<-sum(voters)
    
    # m1 :
    winner_m1<-ifelse(py_voters>total_voters/2, "Py", "R")
    prop_happy_people_m1<-(ifelse(py_voters>total_voters/2, py_voters, total_voters - py_voters)/total_voters)*100
    
    #m2 :
    states_winners <- data.frame(
      s1=ifelse(voters[1]>n/2, "Py", "R"),
      s2=ifelse(voters[2]>n/2, "Py", "R"),
      s3=ifelse(voters[3]>n/2, "Py", "R"))
    winner_m2 <- names(which.max(table(factor(unlist(states_winners), levels = c("Py", "R")))))
    prop_happy_people_m2<-(ifelse(winner_m2=="Py", py_voters, total_voters-py_voters)/total_voters)*100
    
    prop_happy_voters[b,]<-c(prop_happy_people_m1, prop_happy_people_m2)
  }
  
  final_results <- colMeans(prop_happy_voters)
  
  return(final_results)
}


simul_elec(n=101101, case="IC", B=1000)

####
