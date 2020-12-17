require(tidyverse)
require(patchwork)
require(Cairo)

data_wd <- "Simulation_Results"
plot_wd <- "Simulation_Plots"

load(paste0(data_wd,"/ResFull.Rda"))

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
NumSims <- 250

#### Analysis and Plotting of Simulation Results (Both Full Color and AJE versions):
for (k.par in disp_params) {
  ResFull_k <- ResFull %>% filter(k==k.par)
  if (k.par==disp_params[1]) {
    FigNo <- c("1","2")
  } else if (k.par==disp_params[2]) {
    FigNo <- c("S1","S2")
  } else {
    FigNo <- c("S3","S4")
  }

  ### Estimation Figures:
  Infections <- ResFull_k %>% group_by(Time,Correlation) %>% summarize(NSims=n(),
                                                                   median_S1=median(S1Infs, na.rm=TRUE), lower_S1=quantile(S1Infs, 0.25, na.rm=TRUE), upper_S1=quantile(S1Infs, 0.75, na.rm=TRUE),
                                                                   median_S2=median(S2Infs, na.rm=TRUE), lower_S2=quantile(S2Infs, 0.25, na.rm=TRUE), upper_S2=quantile(S2Infs, 0.75, na.rm=TRUE), 
                                                                   median_S3=median(S3Infs, na.rm=TRUE), lower_S3=quantile(S3Infs, 0.25, na.rm=TRUE), upper_S3=quantile(S3Infs, 0.75, na.rm=TRUE))
  Infections2 <- Infections %>% pivot_longer(cols=starts_with(c("median_","upper_","lower_")), names_to=c("value_type","Method"), values_to="values", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")
  Infections2$X <- Infections2$Time + ifelse(Infections2$Method=="S1",-1,ifelse(Infections2$Method=="S2",0,ifelse(Infections2$Method=="S3",1,0)))
  g_inf <- ggplot(Infections2, aes(x=X, y=median, color=Method, shape=Method)) + 
    geom_point(size=2) + geom_errorbar(aes(ymin=lower, ymax=upper)) +
    facet_wrap(facets=vars(Correlation), nrow=2, ncol=2, scales="fixed") +
    scale_x_continuous(name="Time", breaks=samp_times, labels=samp_times) +
    scale_y_continuous(name="Infections", breaks=seq(0,400,by=100), labels=seq(0,400,by=100)) +
    coord_cartesian(ylim=c(0,475), clip="off") +
    scale_color_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error")) + 
    scale_shape_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error"))
  g_inf
  
  AJE_base <- tibble(label=rep(c("Day 25","  Regular Sampling","  Snowball Sampling","  Snowball, Contact Error",
                                 "Day 35","  Regular Sampling","  Snowball Sampling","  Snowball, Contact Error",
                                 "Day 45","  Regular Sampling","  Snowball Sampling","  Snowball, Contact Error"),3),
                     Y=rep(c(10,9.5,8.5,7.5,6.5,6,5,4,3,2.5,1.5,0.5),3), Correlation=rep(corrs, each=12),
                     Time=rep(c(NA,25,25,25,NA,35,35,35,NA,45,45,45), 3),
                     Method=rep(c(NA,"S1","S2","S3"),9))
  AJE_base$ICC <- paste0("ICC=",format(AJE_base$Correlation, nsmall=2))
  Infections2AJE <- full_join(AJE_base, Infections2 %>% dplyr::select(Time,Correlation,Method,median,upper,lower))
  Infections2AJE$IQR <- Infections2AJE$upper - Infections2AJE$lower
  Infections2AJE$Vals <- ifelse(is.na(Infections2AJE$median),NA,
                                paste0(format(Infections2AJE$median,digits=0)," (",
                                       format(Infections2AJE$lower,digits=0),"\U2013",
                                       format(Infections2AJE$upper,digits=0),")"))
  
  p.1.a <- ggplot(data=Infections2AJE, aes(x=median, y=Y)) +
    geom_point(size=2.5) + geom_errorbarh(aes(xmin=lower, xmax=upper, y=Y), inherit.aes=FALSE) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=11, color="black"),
          axis.line.x=element_blank(),
          strip.text=element_text(size=11),
          axis.ticks.x=element_line(color="black", size=1), 
          text=element_text(size=11, color="black"),
          plot.tag=element_text(size=11, color="black"),
          plot.title=element_blank()) +
    labs(tag="A)") +
    facet_wrap(~ICC, nrow=1, ncol=3, scales="free_x") +
    scale_x_continuous(name="Infections", limit=c(-500,800), ## Limit is big enough to allow the text labels
                       breaks=seq(0,500,by=100), labels=as.character(seq(0,500,by=100))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,10.8), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-500, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-250, y=10.4, label='underline("Sampling Approach")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) + 
    geom_text(aes(x=500, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=650, y=10.4, label='underline("Median (IQR)")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) +
    annotate("segment", x=-3, xend=503, y=0, yend=0, color="black", size=2) ## Bottom axis that ends at break labels
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[1],"_A.eps"),
         plot=p.1.a, width=15, height=6, units="in", device=cairo_ps)
  
  Prop.Est <- ResFull_k %>% group_by(Time, Correlation) %>% summarize(NSims=n(),
                                                                    mean_T=mean(TrueSympPerc, na.rm=TRUE), lower_T=mean(TrueSympPerc, na.rm=TRUE)-sd(TrueSympPerc, na.rm=TRUE), upper_T=mean(TrueSympPerc, na.rm=TRUE)+sd(TrueSympPerc, na.rm=TRUE),
                                                                    mean_S1=mean(S1.prop.Est, na.rm=TRUE), lower_S1=mean(S1.prop.Est, na.rm=TRUE)-sd(S1.prop.Est, na.rm=TRUE), upper_S1=mean(S1.prop.Est, na.rm=TRUE)+sd(S1.prop.Est, na.rm=TRUE),
                                                                    mean_S2=mean(S2.prop.Est, na.rm=TRUE), lower_S2=mean(S2.prop.Est, na.rm=TRUE)-sd(S2.prop.Est, na.rm=TRUE), upper_S2=mean(S2.prop.Est, na.rm=TRUE)+sd(S2.prop.Est, na.rm=TRUE),
                                                                    mean_S3=mean(S3.prop.Est, na.rm=TRUE), lower_S3=mean(S3.prop.Est, na.rm=TRUE)-sd(S3.prop.Est, na.rm=TRUE), upper_S3=mean(S3.prop.Est, na.rm=TRUE)+sd(S3.prop.Est, na.rm=TRUE))
  Prop.Est.quarts <- ResFull_k %>% group_by(Time, Correlation) %>% summarize(NSims=n(),
                                                                           median_T=median(TrueSympPerc, na.rm=TRUE), lower_T=quantile(TrueSympPerc, 0.25, na.rm=TRUE), upper_T=quantile(TrueSympPerc, 0.75, na.rm=TRUE),
                                                                           median_S1=median(S1.prop.Est, na.rm=TRUE), lower_S1=quantile(S1.prop.Est, 0.25, na.rm=TRUE), upper_S1=quantile(S1.prop.Est, 0.75, na.rm=TRUE),
                                                                           median_S2=median(S2.prop.Est, na.rm=TRUE), lower_S2=quantile(S2.prop.Est, 0.25, na.rm=TRUE), upper_S2=quantile(S2.prop.Est, 0.75, na.rm=TRUE),
                                                                           median_S3=median(S3.prop.Est, na.rm=TRUE), lower_S3=quantile(S3.prop.Est, 0.25, na.rm=TRUE), upper_S3=quantile(S3.prop.Est, 0.75, na.rm=TRUE))
  Prop.Est2 <- Prop.Est %>% pivot_longer(cols=starts_with(c("mean_","upper_","lower_")), names_to=c("value_type","Method"), values_to="values", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")
  Prop.Est.q2 <- Prop.Est.quarts %>% pivot_longer(cols=starts_with(c("median_","upper_","lower_")), names_to=c("value_type","Method"), values_to="values", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")
  Prop.Est.q2$X <- Prop.Est.q2$Time + ifelse(Prop.Est.q2$Method=="S1",-1,ifelse(Prop.Est.q2$Method=="S2",0,ifelse(Prop.Est.q2$Method=="S3",1,2)))
  
  g_prop_est <- ggplot(Prop.Est.q2, aes(x=X, y=median*100, color=Method, shape=Method)) + 
    geom_point(size=2) + geom_errorbar(aes(ymin=lower*100, ymax=upper*100)) +
    geom_hline(yintercept=100*(symp_probs["High"]*high_risk_prob + symp_probs["Low"]*(1-high_risk_prob)), color="red", linetype="dotted") +
    facet_wrap(facets=vars(Correlation), nrow=2, ncol=2, scales="fixed") +
    scale_x_continuous(name="Time", breaks=samp_times, labels=samp_times) +
    scale_y_continuous(name="Estimated Symptom Rate (%)") +
    scale_color_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth")) + 
    scale_shape_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth"))
  g_prop_est
  
  Prop.Est.AJE <- full_join(AJE_base, Prop.Est.q2 %>% dplyr::select(Time,Correlation,Method,median,upper,lower) %>% filter(Method != "T"))
  Prop.Est.AJE$IQR <- Prop.Est.AJE$upper - Prop.Est.AJE$lower
  Prop.Est.AJE$Vals <- ifelse(is.na(Prop.Est.AJE$median),NA,
                                paste0(format(round(Prop.Est.AJE$median,3),nsmall=3)," (",
                                       format(round(Prop.Est.AJE$lower,3),nsmall=3),"\U2013",
                                       format(round(Prop.Est.AJE$upper,3),nsmall=3),")"))
  
  p.1.b <- ggplot(data=Prop.Est.AJE, aes(x=median, y=Y)) +
    geom_point(size=2.5) + geom_errorbarh(aes(xmin=lower, xmax=upper, y=Y), inherit.aes=FALSE) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=11, color="black"),
          axis.line.x=element_blank(),
          strip.text=element_text(size=11),
          axis.ticks.x=element_line(color="black", size=1), 
          text=element_text(size=11, color="black"),
          plot.tag=element_text(size=11, color="black"),
          plot.title=element_blank()) +
    labs(tag="B)") +
    facet_wrap(~ICC, nrow=1, ncol=3, scales="free_x") +
    scale_x_continuous(name="Estimated Symptom Rate", limit=c(-.1,.16), ## Limit is big enough to allow the text labels
                       breaks=seq(0,.08,by=.02), labels=as.character(seq(0,.08,by=.02))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,10.8), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-.1, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-.05, y=10.4, label='underline("Sampling Approach")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) + 
    geom_text(aes(x=.081, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=.12, y=10.4, label='underline("Median (IQR)")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) +
    annotate("segment", x=-.0005, xend=.0805, y=0, yend=0, color="black", size=2) + ## Bottom axis that ends at break labels
    geom_vline(xintercept=.05)
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[1],"_B.eps"),
         plot=p.1.b, width=15, height=6, units="in", device=cairo_ps)
  
  Prop.MSE <- ResFull_k %>% group_by(Time, Correlation) %>% summarize(NSims=n(),
                                                                    MSE_S1=mean(S1.prop.Err^2, na.rm=TRUE),
                                                                    MSE_S2=mean(S2.prop.Err^2, na.rm=TRUE),
                                                                    MSE_S3=mean(S3.prop.Err^2, na.rm=TRUE))
  Prop.MSE2 <- Prop.MSE %>% pivot_longer(cols=starts_with("MSE_"), names_to=c("value_type","Method"), values_to="MSE", names_sep="_")
  Prop.MSE2$RMSE_perc <- sqrt(Prop.MSE2$MSE)*100
  Prop.MSE2$RMSE <- Prop.MSE2$RMSE_perc/100
  Prop.MSE2$X <- Prop.MSE2$Time + ifelse(Prop.MSE2$Method=="S1",-1,ifelse(Prop.MSE2$Method=="S2",0,ifelse(Prop.MSE2$Method=="S3",1,0)))
  
  g_prop_MSE <- ggplot(Prop.MSE2, aes(x=X, y=RMSE_perc, color=Method, shape=Method)) +
    geom_point(size=2) +
    facet_wrap(facets=vars(Correlation), nrow=2, ncol=2, scales="fixed") +
    scale_x_continuous(name="Time", breaks=samp_times, labels=samp_times) +
    scale_y_continuous(name="RMSE (Percentage Points)", breaks=seq(0,6,by=2)) +
    coord_cartesian(ylim=c(0,6.5), clip="off") +
    scale_color_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth")) + 
    scale_shape_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth"))
  g_prop_MSE
  
  Prop.MSE.AJE <- full_join(AJE_base, Prop.MSE2 %>% dplyr::select(Time,Correlation,Method,RMSE))
  Prop.MSE.AJE$Vals <- ifelse(is.na(Prop.MSE.AJE$RMSE),NA,
                              paste0(format(round(Prop.MSE.AJE$RMSE,4),nsmall=4)))
  
  p.1.c <- ggplot(data=Prop.MSE.AJE, aes(x=RMSE, y=Y)) +
    geom_point(size=2.5) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=11, color="black"),
          axis.line.x=element_blank(),
          strip.text=element_text(size=11),
          axis.ticks.x=element_line(color="black", size=1), 
          text=element_text(size=11, color="black"),
          plot.tag=element_text(size=11, color="black"),
          plot.title=element_blank()) +
    labs(tag="C)") +
    facet_wrap(~ICC, nrow=1, ncol=3, scales="free_x") +
    scale_x_continuous(name="Root Mean Squared Error of Estimated Symptom Rate", limit=c(-.05,.06), ## Limit is big enough to allow the text labels
                       breaks=seq(0,.05,by=.01), labels=as.character(seq(0,.05,by=.01))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,10.8), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-.05, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-.025, y=10.4, label='underline("Sampling Approach")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) + 
    geom_text(aes(x=.05, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=.055, y=10.4, label='underline("Value")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) +
    annotate("segment", x=-.0005, xend=.0505, y=0, yend=0, color="black", size=2) ## Bottom axis that ends at break labels
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[1],"_C.eps"),
         plot=p.1.c, width=15, height=6, units="in", device=cairo_ps)
  
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[1],".png"),
         plot=p.1.a+p.1.b+p.1.c+plot_layout(nrow=3, ncol=1), 
         width=15, height=18, units="in", dpi=300)
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[1],".eps"),
         plot=p.1.a+p.1.b+p.1.c+plot_layout(nrow=3, ncol=1), 
         width=15, height=18, units="in", device=cairo_ps)
  
  
  ### Precision Figures:
  Prop.CIW <- ResFull_k %>% group_by(Time, Correlation) %>% summarize(NSims=n(),
                                                                      median_S1=median(S1.prop.CIW, na.rm=TRUE), lower_S1=quantile(S1.prop.CIW, 0.25, na.rm=TRUE), upper_S1=quantile(S1.prop.CIW, 0.75, na.rm=TRUE),
                                                                      median_S2=median(S2.prop.CIW, na.rm=TRUE), lower_S2=quantile(S2.prop.CIW, 0.25, na.rm=TRUE), upper_S2=quantile(S2.prop.CIW, 0.75, na.rm=TRUE),
                                                                      median_S3=median(S3.prop.CIW, na.rm=TRUE), lower_S3=quantile(S3.prop.CIW, 0.25, na.rm=TRUE), upper_S3=quantile(S3.prop.CIW, 0.75, na.rm=TRUE))
  Prop.CIW2 <- Prop.CIW %>% pivot_longer(cols=starts_with(c("median_","upper_","lower_")), names_to=c("value_type","Method"), values_to="values", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")
  Prop.CIW2$X <- Prop.CIW2$Time + ifelse(Prop.CIW2$Method=="S1",-1,ifelse(Prop.CIW2$Method=="S2",0,ifelse(Prop.CIW2$Method=="S3",1,0)))
  
  g_prop_CIW <- ggplot(Prop.CIW2, aes(x=X, y=median, color=Method, shape=Method)) +
    geom_point(size=2) + geom_errorbar(aes(ymin=lower, ymax=upper)) +
    facet_wrap(facets=vars(Correlation), nrow=2, ncol=2, scales="fixed") +
    scale_x_continuous(name="Time", breaks=samp_times, labels=samp_times) +
    scale_y_continuous(name="95% CI Width", breaks=seq(0,.4,by=.1)) +
    coord_cartesian(ylim=c(0,.4), clip="on") +
    scale_color_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth")) + 
    scale_shape_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth"))
  g_prop_CIW
  
  Prop.CIW.AJE <- full_join(AJE_base, Prop.CIW2 %>% dplyr::select(Time,Correlation,Method,median,upper,lower))
  Prop.CIW.AJE$Vals <- ifelse(is.na(Prop.CIW.AJE$median),NA,
                              paste0(format(round(Prop.CIW.AJE$median,3),nsmall=3)," (",
                                     format(round(Prop.CIW.AJE$lower,3),nsmall=3),"\U2013",
                                     format(round(Prop.CIW.AJE$upper,3),nsmall=3),")"))
  Prop.CIW.AJE$star <- ifelse(Prop.CIW.AJE$upper > 0.5,1,0)
  
  p.2.a <- ggplot(data=Prop.CIW.AJE, aes(x=median, y=Y)) +
    geom_point(size=2.5) + geom_errorbarh(aes(xmin=lower, xmax=pmin(upper,0.5), y=Y), inherit.aes=FALSE) +
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=11, color="black"),
          axis.line.x=element_blank(),
          strip.text=element_text(size=11),
          axis.ticks.x=element_line(color="black", size=1), 
          text=element_text(size=11, color="black"),
          plot.tag=element_text(size=11, color="black"),
          plot.title=element_blank()) +
    labs(tag="A)") +
    facet_wrap(~ICC, nrow=1, ncol=3, scales="free_x") +
    scale_x_continuous(name="95% CI Width", limit=c(-0.7,1), ## Limit is big enough to allow the text labels
                       breaks=seq(0,0.5,by=.1), labels=as.character(seq(0,0.5,by=.1))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,10.8), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-0.7, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-.35, y=10.4, label='underline("Sampling Approach")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) + 
    geom_text(aes(x=.505, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=.75, y=10.4, label='underline("Median (IQR)")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) +
    annotate("segment", x=-.005, xend=.5005, y=0, yend=0, color="black", size=2) ## Bottom axis that ends at break labels
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[2],"_A.eps"),
         plot=p.2.a, width=15, height=6, units="in", device=cairo_ps)
  
  
  Prop.CICov <- ResFull_k %>% group_by(Time, Correlation) %>% summarize(NSims=n(),
                                                                      CovIn_S1=mean(S1.prop.cvgIn, na.rm=TRUE),
                                                                      CovIn_S2=mean(S2.prop.cvgIn, na.rm=TRUE),
                                                                      CovIn_S3=mean(S3.prop.cvgIn, na.rm=TRUE),
                                                                      CovOut_S1=mean(S1.prop.cvgOut, na.rm=TRUE),
                                                                      CovOut_S2=mean(S2.prop.cvgOut, na.rm=TRUE),
                                                                      CovOut_S3=mean(S3.prop.cvgOut, na.rm=TRUE))
  Prop.CICov2 <- Prop.CICov %>% pivot_longer(cols=starts_with("Cov"), names_to=c("value_type","Method"), values_to="Coverage", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="Coverage")
  Prop.CICov2$X <- Prop.CICov2$Time + ifelse(Prop.CICov2$Method=="S1",-1,ifelse(Prop.CICov2$Method=="S2",0,ifelse(Prop.CICov2$Method=="S3",1,0)))
  
  g_prop_CICov <- ggplot(Prop.CICov2, aes(x=X, y=CovIn*100, color=Method, shape=Method)) +
    geom_point(size=2) +
    geom_hline(yintercept=95, color="red", linetype="dotted") +
    facet_wrap(facets=vars(Correlation), nrow=2, ncol=2, scales="fixed") +
    scale_x_continuous(name="Time", breaks=samp_times, labels=samp_times) +
    scale_y_continuous(name="95% CI Coverage (%)", breaks=seq(75,100,by=5)) +
    coord_cartesian(ylim=c(75,100), clip="off") +
    scale_color_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth")) + 
    scale_shape_discrete(labels=c("Regular Sampling","Snowball Sampling","Snowball Sampling with Contact Error","Truth"))
  g_prop_CICov
  
  Prop.Cov.AJE <- full_join(AJE_base, Prop.CICov2 %>% dplyr::select(Time,Correlation,Method,CovIn))
  Prop.Cov.AJE$Cov <- Prop.Cov.AJE$CovIn*100
  Prop.Cov.AJE$Vals <- ifelse(is.na(Prop.Cov.AJE$Cov),NA,
                              paste0(format(round(Prop.Cov.AJE$Cov,1),nsmall=1),"%"))
  
  p.2.b <- ggplot(data=Prop.Cov.AJE, aes(x=Cov, y=Y)) +
    geom_point(size=2.5) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=11, color="black"),
          axis.line.x=element_blank(),
          axis.ticks.x=element_line(color="black", size=1), 
          text=element_text(size=11, color="black"),
          strip.text=element_text(size=11),
          plot.tag=element_text(size=11, color="black"),
          plot.title=element_blank()) +
    labs(tag="B)") +
    facet_wrap(~ICC, nrow=1, ncol=3, scales="free_x") +
    scale_x_continuous(name="95% CI Coverage (%)", limit=c(63,105), ## Limit is big enough to allow the text labels
                       breaks=seq(80,100,by=5), labels=as.character(seq(80,100,by=5))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,10.8), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=63, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=71.5, y=10.4, label='underline("Sampling Approach")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) + 
    geom_text(aes(x=101, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=11/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=102.5, y=10.4, label='underline("Coverage")', parse=TRUE, hjust="center", size=11/ggplot2:::.pt) +
    annotate("segment", x=79.95, xend=100.05, y=0, yend=0, color="black", size=2) + ## Bottom axis that ends at break labels
    geom_vline(xintercept=95)
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[2],"_B.eps"),
         plot=p.2.b, width=15, height=6, units="in", device=cairo_ps)
  
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[2],".png"),
         plot=p.2.a+p.2.b+plot_layout(nrow=2, ncol=1), 
         width=15, height=12, units="in", dpi=300)
  
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[2],".eps"),
         plot=p.2.a+p.2.b+plot_layout(nrow=2, ncol=1), 
         width=15, height=12, units="in", device=cairo_ps)
  
  
  ggsave(filename=paste0(plot_wd,"/Full Color/Simulation_Ests_",k.par,".png"),
         plot=g_inf+theme(legend.position="none") + labs(title="A. Number of Infections in Sample") + 
           guide_area() + 
           g_prop_est + theme(legend.position="right") + labs(title="B. Median, Quartiles of Est. Symptom Rate") +
           g_prop_MSE + theme(legend.position="none") + labs(title="C. RMSE of Est. Symptom Rate") + 
           plot_layout(nrow=2, ncol=2, byrow=TRUE, guides="collect"),
         width=8, height=8, units="in", dpi=300)
  
  ggsave(filename=paste0(plot_wd,"/Full Color/Simulation_CIs_",k.par,".png"),
         plot=g_prop_CIW+theme(legend.position="bottom")+labs(title="A. Mean 95% CI Width") +
           g_prop_CICov+theme(legend.position="none")+labs(title="B. Empirical 95% CI Coverage") +
           plot_layout(nrow=1, ncol=2, byrow=TRUE, guides="collect") & theme(legend.position="bottom"),
         width=8, height=5, units="in", dpi=300)
}

## Proportion of Simulations with at least one symptomatic case, by analysis method
Prop.Obs <- ResFull %>% group_by(Time,Correlation,k) %>% 
  summarize(Prop.S1=mean(S1Symps>0), Prop.S2=mean(S2Symps>0), Prop.S3=mean(S3Symps>0))
Prop.Obs.t <- ResFull %>% group_by(Time) %>% 
  summarize(Prop.S1=mean(S1Symps>0), Prop.S2=mean(S2Symps>0), Prop.S3=mean(S3Symps>0))


