require(igraph)
require(tidyverse)
require(survey)
require(patchwork)
require(Cairo)

setwd("Simulation_Results")

## expit function
expit <- function(x) {
  return(1/(1+exp(-x)))
}

## Population Parameters
pop_n <- 10000
avg_contacts <- 20
high_risk_prob <- 0.20

## Infection Parameters
R0 <- 2.5 #Mean of transmission
disp_params <- c(.1,.5,1) #Dispersion parameter for transmission
infday_opts <- 2:9
infday_prob <- rep(1/8,8) #If infected from an infector, 1/8 chance of being infected each day after source's infection, starting 2 days later and ending 9 days later
i0 <- 6 #Number of people infected on day 1

## Symptomatic Parameters
symp_probs <- c("High"=.016/high_risk_prob, "Low"=.034/(1-high_risk_prob))
overall_symp <- symp_probs["High"]*high_risk_prob+symp_probs["Low"]*(1-high_risk_prob)
OR_symp <- (symp_probs["High"]/(1-symp_probs["High"]))/(symp_probs["Low"]/(1-symp_probs["Low"]))
corrs <- c(0,0.05,0.1)

## Sampling Parameters
samp_times <- c(25,35,45) #Last day that will seroconvert
num_samps <- 600 #Number of samples for SRS
num_indices <- 30 #Number of index individuals for snowball sampling
missed_contacts <- 2 #For error in snowball sampling, number of contacts missed by each index
false_contacts <- 2 #For error in snowball sampling, number of extra contacts named by each index


## Set up population
st <- proc.time()
set.seed(19730105)
NumSims <- 250
ResFull <- NULL

for (k.par in disp_params) {
  print(paste0("Starting Dispersion Parameter: ",k.par))
  p.param <- avg_contacts/pop_n
  q.param <- R0/avg_contacts
  var.q <- (1/k.par+1/pop_n)*q.param^2/(p.param*(pop_n-1))
  
  nu.param <- q.param*(1-q.param)/var.q-1
  alpha.param <- nu.param*q.param
  beta.param <- nu.param*(1-q.param)
  
  for (corr in corrs) {
    print(paste0("Starting Correlation: ",corr,", within k=",k.par))
    var.c <- corr*symp_probs*(1-symp_probs)
    nu.c <- symp_probs*(1-symp_probs)/var.c-1
    alpha.c <- symp_probs*nu.c
    beta.c <- nu.c-alpha.c
    for (SimNo in 1:NumSims) {
      st1 <- proc.time()
      print(paste0("Starting Sim Number: ",SimNo,", within k=",k.par,", corr=",corr))
      g <- erdos.renyi.game(n=pop_n,
                            p.or.m=avg_contacts/(pop_n-1),
                            type="gnp",
                            directed=FALSE)

      V(g)$risk <- ifelse(rbinom(n=pop_n, size=1, prob=high_risk_prob)==1,"High","Low")
      V(g)$infday <- 100000
      V(g)$infector <- NA
      V(g)$symp <- NA
      V(g)$contacts <- degree(g)

      ## Start infections
      StartInfs <- as.numeric(sample(V(g), size=i0, replace=FALSE))
      V(g)$infday[StartInfs] <- 1
      V(g)$infector[StartInfs] <- 0
      V(g)$symp[StartInfs] <- rbinom(n=1, size=1, prob=symp_probs[V(g)$risk[StartInfs]])

      ## Time steps
      times <- 1:max(samp_times)
      for (t in times) {
        inf_nodes <- as.numeric(V(g)[V(g)$infday==t])
        contacts <- V(g)[inf_nodes]$contacts
        qs <- rbeta(n=length(inf_nodes), shape1=alpha.param, shape2=beta.param)
        infs <- rbinom(n=length(inf_nodes), size=contacts, prob=qs)
        new_infs <- NULL
        for (j in seq_along(inf_nodes)) {
          if (infs[j] > 0) {
            if (corr==0) {
              symp_probs_j <- symp_probs
            } else {
              symp_probs_j <- c("High"=rbeta(1, alpha.c["High"], beta.c["High"]),
                                "Low"=rbeta(1, alpha.c["Low"], beta.c["Low"]))
            }
            new_infs <- sample(as.numeric(adjacent_vertices(g, inf_nodes[j])[[1]]), infs[j], replace=FALSE)
            new_days <- sample(x=infday_opts, size=length(new_infs), replace=TRUE, prob=infday_prob)
            for (k in seq_along(new_infs)) {
              if (t+new_days[k] < V(g)$infday[new_infs[k]]) {
                V(g)$infday[new_infs[k]] <- t+new_days[k]
                V(g)$infector[new_infs[k]] <- inf_nodes[j]
                V(g)$symp[new_infs[k]] <- rbinom(n=1, size=1, prob=symp_probs_j[V(g)$risk[StartInfs]])
              }
              else if (t+new_days[k] == V(g)$infday[new_infs[k]]) {
                if (rbinom(1, 1, .5)==1) {
                  V(g)$infday[new_infs[k]] <- t+new_days[k]
                  V(g)$infector[new_infs[k]] <- inf_nodes[j]
                  V(g)$symp[new_infs[k]] <- rbinom(n=1, size=1, prob=symp_probs_j[V(g)$risk[StartInfs]])
                }
              }
            }
          }
        }
      }
      
      for (samp_time in samp_times) {
        ## Summarizing Results of Outbreak
        Res <- igraph::as_data_frame(g, what="vertices")
        Res$infected <- ifelse(Res$infday <= max(times), 1, 0)
        Res <- as_tibble(Res)
        Res$vertex <- as.numeric(V(g))
        
        Res_InfDays <- Res %>% filter(infected==1) %>% group_by(infday) %>% 
          summarize(num_inf=n(), num_symp=sum(symp))
        Res_InfDays$CI <- cumsum(Res_InfDays$num_inf)
        Res_InfDays$SympCI <- cumsum(Res_InfDays$num_symp)
        Res_SympRisk <- Res %>% filter(infday <= samp_time) %>% group_by(risk,symp) %>% 
          summarize(num=n())
        
        ## Regular Sampling
        samp1 <- slice_sample(Res, n=num_samps, replace=FALSE)
        samp1$infected <- ifelse(samp1$infday <= samp_time, 1, 0)
        samp1I <- samp1 %>% filter(infected==1)
        if (dim(samp1I)[1] > 0) {
          Reg1a <- summary(glm(symp~1, family=binomial, data=samp1I))
          Reg1aRes <- expit(c(Reg1a$coefficients["(Intercept)","Estimate"]-1.96*Reg1a$coefficients["(Intercept)","Std. Error"],
                            Reg1a$coefficients["(Intercept)","Estimate"],
                            Reg1a$coefficients["(Intercept)","Estimate"]+1.96*Reg1a$coefficients["(Intercept)","Std. Error"]))
          if (sum(samp1I$risk=="High") > 0 & sum(samp1I$risk=="Low") > 0) {
            Reg1b <- summary(glm(symp~risk, family=binomial, data=samp1I))
            Reg1bRes <- pmin(exp(c(-Reg1b$coefficients["riskLow","Estimate"]-1.96*Reg1b$coefficients["riskLow","Std. Error"],
                              -Reg1b$coefficients["riskLow","Estimate"],
                              -Reg1b$coefficients["riskLow","Estimate"]+1.96*Reg1b$coefficients["riskLow","Std. Error"])),100000)
          } else {
            Reg1bRes <- c(NA,NA,NA)
          }
        } else {
          Reg1aRes <- c(NA,NA,NA)
          Reg1bRes <- Reg1aRes
        }
  
        ## Snowball Sampling
        samp2ind <- slice_sample(Res %>% filter(infday <= samp_time-max(infday_opts)), 
                                 n=num_indices, replace=FALSE)
        samp2 <- NULL
        for (i in 1:length(samp2ind$vertex)) {
          samp2 <- bind_rows(samp2, bind_cols(Index=samp2ind$vertex[i],
                                              vertex=unlist(adjacent_vertices(g, samp2ind$vertex[i]))))
        }
        samp2 <- samp2 %>% distinct(vertex, .keep_all=TRUE) %>% filter(!(vertex %in% samp2ind$vertex)) ## Keeps the first appearance of each vertex (by index individual, where index order is random)
        samp2 <- left_join(samp2, Res, by="vertex")
        samp2$infected <- ifelse(samp2$infday <= samp_time, 1, 0)
        samp2I <- samp2 %>% filter(infected==1)
        if (dim(samp2I)[1] > 0) {
          svydes <- svydesign(id=~Index, data=samp2I)
          Reg2a <- summary(svyglm(symp~1, design=svydes, family=binomial))
          Reg2aRes <- expit(c(Reg2a$coefficients["(Intercept)","Estimate"]-1.96*Reg2a$coefficients["(Intercept)","Std. Error"],
                            Reg2a$coefficients["(Intercept)","Estimate"],
                            Reg2a$coefficients["(Intercept)","Estimate"]+1.96*Reg2a$coefficients["(Intercept)","Std. Error"]))
          if (sum(samp2I$risk=="High") > 0 & sum(samp2I$risk=="Low") > 0) {
            Reg2b <- summary(svyglm(symp~risk, design=svydes, family=binomial))
            Reg2bRes <- pmin(exp(c(-Reg2b$coefficients["riskLow","Estimate"]-1.96*Reg2b$coefficients["riskLow","Std. Error"],
                                   -Reg2b$coefficients["riskLow","Estimate"],
                                   -Reg2b$coefficients["riskLow","Estimate"]+1.96*Reg2b$coefficients["riskLow","Std. Error"])),100000)
          } else {
            Reg2bRes <- c(NA,NA,NA)
          }
        } else {
          Reg2aRes <- c(NA,NA,NA)
          Reg2bRes <- Reg2aRes
        }
        ##Snowball Sampling with Errors in Contacts
        samp3 <- NULL
        for (i in 1:length(samp2ind$vertex)) {
          verts <- unlist(adjacent_vertices(g, samp2ind$vertex[i]))
          verts.add <- sample(setdiff(as.numeric(V(g)),verts), false_contacts)
          verts.adj <- c(sample(verts, length(verts)-missed_contacts),verts.add)
          samp3 <- bind_rows(samp3, bind_cols(Index=samp2ind$vertex[i],
                                              vertex=verts.adj))
        }
        samp3 <- samp3 %>% distinct(vertex, .keep_all=TRUE) %>% filter(!(vertex %in% samp2ind$vertex)) ## Keeps the first appearance of each vertex (by index individual, where index order is random)
        samp3 <- left_join(samp3, Res, by="vertex")
        samp3$infected <- ifelse(samp3$infday <= samp_time, 1, 0)
        samp3I <- samp3 %>% filter(infected==1)
        if (dim(samp3I)[1] > 0) {
          svydes3 <- svydesign(id=~Index, data=samp3I)
          Reg3a <- summary(svyglm(symp~1, design=svydes3, family=binomial))
          Reg3aRes <- expit(c(Reg3a$coefficients["(Intercept)","Estimate"]-1.96*Reg3a$coefficients["(Intercept)","Std. Error"],
                            Reg3a$coefficients["(Intercept)","Estimate"],
                            Reg3a$coefficients["(Intercept)","Estimate"]+1.96*Reg3a$coefficients["(Intercept)","Std. Error"]))
          if (sum(samp3I$risk=="High") > 0 & sum(samp3I$risk=="Low") > 0) {
            Reg3b <- summary(svyglm(symp~risk, design=svydes3, family=binomial))
            Reg3bRes <- pmin(exp(c(-Reg3b$coefficients["riskLow","Estimate"]-1.96*Reg3b$coefficients["riskLow","Std. Error"],
                                   -Reg3b$coefficients["riskLow","Estimate"],
                                   -Reg3b$coefficients["riskLow","Estimate"]+1.96*Reg3b$coefficients["riskLow","Std. Error"])),100000)
          } else {
            Reg3bRes <- c(NA,NA,NA)
          }
        } else {
          Reg3aRes <- c(NA,NA,NA)
          Reg3bRes <- c(NA,NA,NA)
        }
        ResFull <- bind_rows(ResFull, 
                             tibble(Sim=SimNo, Time=samp_time, Correlation=corr, k=k.par,
                                    CI=Res_InfDays$CI[Res_InfDays$infday==samp_time]/pop_n,
                                    S1Infs=sum(samp1$infected), S1Symps=sum(samp1I$symp, na.rm=TRUE),
                                    S1HRInf=sum(ifelse(samp1I$risk=="High",1,0)), S1HRSymp=sum(samp1I$symp*ifelse(samp1I$risk=="High",1,0), na.rm=TRUE),
                                    S2Infs=sum(samp2$infected), S2Symps=sum(samp2I$symp, na.rm=TRUE),
                                    S2HRInf=sum(ifelse(samp2I$risk=="High",1,0)), S2HRSymp=sum(samp2I$symp*ifelse(samp2I$risk=="High",1,0), na.rm=TRUE),
                                    S2N=length(samp2$infected),
                                    S3Infs=sum(samp3$infected), S3Symps=sum(samp3I$symp, na.rm=TRUE),
                                    S3HRInf=sum(ifelse(samp3I$risk=="High",1,0)), S3HRSymp=sum(samp3I$symp*ifelse(samp3I$risk=="High",1,0), na.rm=TRUE),
                                    S3N=length(samp3$infected),
                                    TrueInfPerc=Res_InfDays$CI[Res_InfDays$infday==samp_time]/pop_n,
                                    TrueSympPerc=Res_InfDays$SympCI[Res_InfDays$infday==samp_time]/Res_InfDays$CI[Res_InfDays$infday==samp_time],
                                    TrueOR=ifelse(Res_SympRisk$num[Res_SympRisk$risk=="Low" & Res_SympRisk==1]==0,NA,
                                                  (Res_SympRisk$num[Res_SympRisk$risk=="High" & Res_SympRisk$symp==1]/Res_SympRisk$num[Res_SympRisk$risk=="High" & Res_SympRisk$symp==0])/
                                                    (Res_SympRisk$num[Res_SympRisk$risk=="Low" & Res_SympRisk$symp==1]/Res_SympRisk$num[Res_SympRisk$risk=="Low" & Res_SympRisk$symp==0])),
                                    S1.prop.Est=Reg1aRes[2], S1.prop.CIL=Reg1aRes[1], S1.prop.CIU=Reg1aRes[3],
                                    S2.prop.Est=Reg2aRes[2], S2.prop.CIL=Reg2aRes[1], S2.prop.CIU=Reg2aRes[3],
                                    S3.prop.Est=Reg3aRes[2], S3.prop.CIL=Reg3aRes[1], S3.prop.CIU=Reg3aRes[3],
                                    S1.OR.Est=Reg1bRes[2], S1.OR.CIL=Reg1bRes[1], S1.OR.CIU=Reg1bRes[3],
                                    S2.OR.Est=Reg2bRes[2], S2.OR.CIL=Reg2bRes[1], S2.OR.CIU=Reg2bRes[3],
                                    S3.OR.Est=Reg3bRes[2], S3.OR.CIL=Reg3bRes[1], S3.OR.CIU=Reg3bRes[3]))
      }
      print(proc.time()-st1)
    }
  }
}
proc.time() - st
save(ResFull, file="ResFull_pre.Rda")

##Post-Processing of Results:
ResFull$S1.prop.CIW <- ResFull$S1.prop.CIU-ResFull$S1.prop.CIL
ResFull$S2.prop.CIW <- ResFull$S2.prop.CIU-ResFull$S2.prop.CIL
ResFull$S3.prop.CIW <- ResFull$S3.prop.CIU-ResFull$S3.prop.CIL
ResFull$S1.OR.CIW <- ResFull$S1.OR.CIU-ResFull$S1.OR.CIL
ResFull$S2.OR.CIW <- ResFull$S2.OR.CIU-ResFull$S2.OR.CIL
ResFull$S3.OR.CIW <- ResFull$S3.OR.CIU-ResFull$S3.OR.CIL
ResFull$S1.prop.cvgIn <- ifelse(ResFull$TrueSympPerc <= ResFull$S1.prop.CIU & ResFull$TrueSympPerc >= ResFull$S1.prop.CIL,1,0)
ResFull$S2.prop.cvgIn <- ifelse(ResFull$TrueSympPerc <= ResFull$S2.prop.CIU & ResFull$TrueSympPerc >= ResFull$S2.prop.CIL,1,0)
ResFull$S3.prop.cvgIn <- ifelse(ResFull$TrueSympPerc <= ResFull$S3.prop.CIU & ResFull$TrueSympPerc >= ResFull$S3.prop.CIL,1,0)
ResFull$S1.prop.cvgOut <- ifelse(overall_symp <= ResFull$S1.prop.CIU & overall_symp >= ResFull$S1.prop.CIL,1,0)
ResFull$S2.prop.cvgOut <- ifelse(overall_symp <= ResFull$S2.prop.CIU & overall_symp >= ResFull$S2.prop.CIL,1,0)
ResFull$S3.prop.cvgOut <- ifelse(overall_symp <= ResFull$S3.prop.CIU & overall_symp >= ResFull$S3.prop.CIL,1,0)
ResFull$S1.OR.cvgIn <- ifelse(is.na(ResFull$TrueOR),NA,
                              ifelse(ResFull$TrueOR <= ResFull$S1.OR.CIU & ResFull$TrueOR >= ResFull$S1.OR.CIL,1,0))
ResFull$S2.OR.cvgIn <- ifelse(is.na(ResFull$TrueOR),NA,
                              ifelse(ResFull$TrueOR <= ResFull$S2.OR.CIU & ResFull$TrueOR >= ResFull$S2.OR.CIL,1,0))
ResFull$S3.OR.cvgIn <- ifelse(is.na(ResFull$TrueOR),NA,
                              ifelse(ResFull$TrueOR <= ResFull$S3.OR.CIU & ResFull$TrueOR >= ResFull$S3.OR.CIL,1,0))
ResFull$S1.OR.cvgOut <- ifelse(OR_symp <= ResFull$S1.OR.CIU & OR_symp >= ResFull$S1.OR.CIL, 1, 0)
ResFull$S2.OR.cvgOut <- ifelse(OR_symp <= ResFull$S2.OR.CIU & OR_symp >= ResFull$S2.OR.CIL, 1, 0)
ResFull$S3.OR.cvgOut <- ifelse(OR_symp <= ResFull$S3.OR.CIU & OR_symp >= ResFull$S3.OR.CIL, 1, 0)
ResFull$S1.prop.Err <- ResFull$S1.prop.Est-ResFull$TrueSympPerc
ResFull$S2.prop.Err <- ResFull$S2.prop.Est-ResFull$TrueSympPerc
ResFull$S3.prop.Err <- ResFull$S3.prop.Est-ResFull$TrueSympPerc

## Save Results:
save(ResFull, file="ResFull.Rda")


