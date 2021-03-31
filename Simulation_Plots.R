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
for (k.par in disp_params) { #disp_params[1]) { 
  ResFull_k <- ResFull %>% filter(k==k.par)
  if (k.par==disp_params[1]) {
    FigNo <- c("1","2","3")
  } else if (k.par==disp_params[2]) {
    FigNo <- c("W1","W2","W3")
  } else {
    FigNo <- c("W4","W5","W6")
  }
  for (icc in corrs) {
    ResFull_use <- ResFull_k %>% filter(Correlation==icc)
    if (icc==corrs[1]) {
      FigLab <- c("A","C","E")
    } else if (icc==corrs[2]) {
      FigLab <- c("B","D","F")
    } else {
      FigLab <- c("A","B","C","D","E")
    }

  ### Estimation Figures:
  Infections <- ResFull_use %>% group_by(Time,Correlation) %>% 
    summarize(NSims=n(),median_S1=median(S1Infs, na.rm=TRUE),
              lower_S1=quantile(S1Infs, 0.25, na.rm=TRUE),upper_S1=quantile(S1Infs, 0.75, na.rm=TRUE),
              median_S2=median(S2Infs, na.rm=TRUE), lower_S2=quantile(S2Infs, 0.25, na.rm=TRUE), 
              upper_S2=quantile(S2Infs, 0.75, na.rm=TRUE), median_S3=median(S3Infs, na.rm=TRUE), 
              lower_S3=quantile(S3Infs, 0.25, na.rm=TRUE), upper_S3=quantile(S3Infs, 0.75, na.rm=TRUE))
  Infections2 <- Infections %>% pivot_longer(cols=starts_with(c("median_","upper_","lower_")), names_to=c("value_type","Method"), values_to="values", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")

  AJE_base <- tibble(label=c("Day 25",paste0("  ","Regular sampling"),
                                 paste0("  ","Snowball sampling"),paste0("  ","Snowball, error"),
                                 "Day 35",paste0("  ","Regular sampling"),
                                 paste0("  ","Snowball sampling"),paste0("  ","Snowball, error"),
                                 "Day 45",paste0("  ","Regular sampling"),
                                 paste0("  ","Snowball sampling"),paste0("  ","Snowball, error")),
                     Y=seq(from=4.8, to=0.4, by=-0.4), Correlation=rep(icc,12),
                     Time=c(NA,25,25,25,NA,35,35,35,NA,45,45,45),
                     Method=rep(c(NA,"S1","S2","S3"),3))
  Infections2AJE <- full_join(AJE_base, Infections2 %>% dplyr::select(Time,Correlation,Method,median,upper,lower))
  Infections2AJE$IQR <- Infections2AJE$upper - Infections2AJE$lower
  Infections2AJE$Vals <- ifelse(is.na(Infections2AJE$median),NA,
                                paste0(format(Infections2AJE$median, digits=0, trim=FALSE)," (",
                                       format(Infections2AJE$lower, digits=0, trim=TRUE),"\U2013",
                                       format(Infections2AJE$upper, digits=0, trim=TRUE),")"))
  
  plot_top <- ggplot(data=Infections2AJE, aes(x=median, y=Y)) +
    geom_point(size=2.5) + geom_errorbarh(aes(xmin=lower, xmax=upper, y=Y), height=.3, inherit.aes=FALSE) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=8, color="black"),
          axis.title.x=element_text(size=8, color="black", hjust=.54),
          axis.line.x=element_blank(),
          strip.text=element_text(size=8),
          axis.ticks.x=element_line(color="black", size=1/ggplot2:::.pt), 
          text=element_text(size=8, color="black"),
          plot.tag=element_text(size=8, color="black"),
          plot.title=element_blank()) +
    labs(tag=paste0(FigLab[1],")")) +
    scale_x_continuous(name="Median", limit=c(-475,850), ## Limit is big enough to allow the text labels
                       breaks=seq(0,500,by=100), labels=as.character(seq(0,500,by=100))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,5.6), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-475, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-245, y=5.2, label='underline("Day and Sampling")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) + 
    geom_text(aes(x=550, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=725, y=5.2, label='underline("Median (IQR)")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) +
    annotate("segment", x=-2.3, xend=502.3, y=0, yend=0, color="black", size=1) ## Bottom axis that ends at break labels
  
  ggsave(filename=paste0(plot_wd,"/Individual_Panels/Fig",ifelse(icc==corrs[3],FigNo[3],FigNo[1]),"_",FigLab[1],".eps"),
         plot=plot_top, width=3.5, height=2.5, units="in", device=cairo_ps)
  
  assign(x=paste0("p.",ifelse(icc==corrs[3],"S","1"),".",FigLab[1]), value=plot_top)
  
  Prop.Est <- ResFull_use %>% group_by(Time,Correlation) %>% 
    summarize(NSims=n(), mean_T=mean(TrueSympPerc, na.rm=TRUE), 
              lower_T=mean(TrueSympPerc, na.rm=TRUE)-sd(TrueSympPerc, na.rm=TRUE), 
              upper_T=mean(TrueSympPerc, na.rm=TRUE)+sd(TrueSympPerc, na.rm=TRUE),
              mean_S1=mean(S1.prop.Est, na.rm=TRUE), 
              lower_S1=mean(S1.prop.Est, na.rm=TRUE)-sd(S1.prop.Est, na.rm=TRUE), 
              upper_S1=mean(S1.prop.Est, na.rm=TRUE)+sd(S1.prop.Est, na.rm=TRUE),
              mean_S2=mean(S2.prop.Est, na.rm=TRUE), 
              lower_S2=mean(S2.prop.Est, na.rm=TRUE)-sd(S2.prop.Est, na.rm=TRUE), 
              upper_S2=mean(S2.prop.Est, na.rm=TRUE)+sd(S2.prop.Est, na.rm=TRUE),
              mean_S3=mean(S3.prop.Est, na.rm=TRUE), 
              lower_S3=mean(S3.prop.Est, na.rm=TRUE)-sd(S3.prop.Est, na.rm=TRUE), 
              upper_S3=mean(S3.prop.Est, na.rm=TRUE)+sd(S3.prop.Est, na.rm=TRUE))
  Prop.Est.quarts <- ResFull_use %>% group_by(Time,Correlation) %>% 
    summarize(NSims=n(), median_T=median(TrueSympPerc, na.rm=TRUE), 
              lower_T=quantile(TrueSympPerc, 0.25, na.rm=TRUE), upper_T=quantile(TrueSympPerc, 0.75, na.rm=TRUE),
              median_S1=median(S1.prop.Est, na.rm=TRUE), lower_S1=quantile(S1.prop.Est, 0.25, na.rm=TRUE), 
              upper_S1=quantile(S1.prop.Est, 0.75, na.rm=TRUE), median_S2=median(S2.prop.Est, na.rm=TRUE), 
              lower_S2=quantile(S2.prop.Est, 0.25, na.rm=TRUE), upper_S2=quantile(S2.prop.Est, 0.75, na.rm=TRUE),
              median_S3=median(S3.prop.Est, na.rm=TRUE), lower_S3=quantile(S3.prop.Est, 0.25, na.rm=TRUE), 
              upper_S3=quantile(S3.prop.Est, 0.75, na.rm=TRUE))
  Prop.Est2 <- Prop.Est %>% pivot_longer(cols=starts_with(c("mean_","upper_","lower_")), 
                                         names_to=c("value_type","Method"), values_to="values", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")
  Prop.Est.q2 <- Prop.Est.quarts %>% pivot_longer(cols=starts_with(c("median_","upper_","lower_")), 
                                                  names_to=c("value_type","Method"), values_to="values", 
                                                  names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")

  Prop.Est.AJE <- full_join(AJE_base, Prop.Est.q2 %>% dplyr::select(Time,Correlation,Method,median,upper,lower) %>% 
                              filter(Method != "T"))
  Prop.Est.AJE$IQR <- Prop.Est.AJE$upper - Prop.Est.AJE$lower
  Prop.Est.AJE$Vals <- ifelse(is.na(Prop.Est.AJE$median),NA,
                                paste0(format(round(Prop.Est.AJE$median,3),nsmall=3, trim=FALSE)," (",
                                       format(round(Prop.Est.AJE$lower,3),nsmall=3, trim=TRUE),"\U2013",
                                       format(round(Prop.Est.AJE$upper,3),nsmall=3, trim=TRUE),")"))
  
  plot_mid <- ggplot(data=Prop.Est.AJE, aes(x=median, y=Y)) +
    geom_point(size=2.5) + geom_errorbarh(aes(xmin=lower, xmax=upper, y=Y), height=.3, inherit.aes=FALSE) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=8, color="black"),
          axis.line.x=element_blank(),
          axis.title.x=element_text(size=8, color="black", hjust=.515),
          strip.text=element_text(size=8),
          axis.ticks.x=element_line(color="black", size=1/ggplot2:::.pt), 
          text=element_text(size=8, color="black"),
          plot.tag=element_text(size=8, color="black"),
          plot.title=element_blank()) +
    labs(tag=paste0(FigLab[2],")")) +
    scale_x_continuous(name="Median", limit=c(-.13,.20), ## Limit is big enough to allow the text labels
                       breaks=seq(0,.08,by=.04), labels=as.character(seq(0,.08,by=.04))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,5.6), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-.13, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-.072, y=5.2, label='underline("Day and Sampling")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) + 
    geom_text(aes(x=.1, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=.16, y=5.2, label='underline("Median (IQR)")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) +
    annotate("segment", x=-.000525, xend=.080525, y=0, yend=0, color="black", size=1) + ## Bottom axis that ends at break labels
    geom_segment(x=.05, xend=.05, y=0, yend=4.8, color="black", size=1/ggplot2:::.pt)
  
  ggsave(filename=paste0(plot_wd,"/Individual_Panels/Fig",ifelse(icc==corrs[3],FigNo[3],FigNo[1]),"_",FigLab[2],".eps"),
         plot=plot_mid, width=3.5, height=2.5, units="in", device=cairo_ps)
  
  assign(x=paste0("p.",ifelse(icc==corrs[3],"S","1"),".",FigLab[2]), value=plot_mid)
  
  Prop.MSE <- ResFull_use %>% group_by(Time, Correlation) %>% 
    summarize(NSims=n(), MSE_S1=mean(S1.prop.Err^2, na.rm=TRUE),
              MSE_S2=mean(S2.prop.Err^2, na.rm=TRUE),  MSE_S3=mean(S3.prop.Err^2, na.rm=TRUE))
  Prop.MSE2 <- Prop.MSE %>% pivot_longer(cols=starts_with("MSE_"), names_to=c("value_type","Method"), 
                                         values_to="MSE", names_sep="_")
  Prop.MSE2$RMSE_perc <- sqrt(Prop.MSE2$MSE)*100
  Prop.MSE2$RMSE <- Prop.MSE2$RMSE_perc/100

  Prop.MSE.AJE <- full_join(AJE_base, Prop.MSE2 %>% dplyr::select(Time,Correlation,Method,RMSE))
  Prop.MSE.AJE$Vals <- ifelse(is.na(Prop.MSE.AJE$RMSE),NA,
                              paste0(format(round(Prop.MSE.AJE$RMSE,4),nsmall=4)))
  
  plot_bot <- ggplot(data=Prop.MSE.AJE, aes(x=RMSE, y=Y)) +
    geom_point(size=2.5) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=8, color="black"),
          axis.line.x=element_blank(),
          axis.title.x=element_text(size=8, color="black", hjust=.67),
          strip.text=element_text(size=8),
          axis.ticks.x=element_line(color="black", size=1/ggplot2:::.pt), 
          text=element_text(size=8, color="black"),
          plot.tag=element_text(size=8, color="black"),
          plot.title=element_blank()) +
    labs(tag=paste0(FigLab[3],")")) +
    scale_x_continuous(name="Root Mean Squared Error", limit=c(-.039,.066), ## Limit is big enough to allow the text labels
                       breaks=seq(0,.05,by=.025), labels=as.character(seq(0,.05,by=.025))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,5.6), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-.039, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-.021, y=5.2, label='underline("Day and Sampling")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) + 
    geom_text(aes(x=.055, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=.055, y=5.2, label='underline("RMSE")', parse=TRUE, hjust="left", size=8/ggplot2:::.pt) +
    annotate("segment", x=-.00018, xend=.05018, y=0, yend=0, color="black", size=1) ## Bottom axis that ends at break labels

  ggsave(filename=paste0(plot_wd,"/Individual_Panels/Fig",ifelse(icc==corrs[3],FigNo[3],FigNo[1]),"_",FigLab[3],".eps"),
         plot=plot_bot, width=3.5, height=2.5, units="in", device=cairo_ps)
  
  assign(x=paste0("p.",ifelse(icc==corrs[3],"S","1"),".",FigLab[3]), value=plot_bot)
  
  ### Precision Figures:
  Prop.CIW <- ResFull_use %>% group_by(Time, Correlation) %>% 
    summarize(NSims=n(), median_S1=median(S1.prop.CIW, na.rm=TRUE), 
              lower_S1=quantile(S1.prop.CIW, 0.25, na.rm=TRUE), upper_S1=quantile(S1.prop.CIW, 0.75, na.rm=TRUE),
              median_S2=median(S2.prop.CIW, na.rm=TRUE), lower_S2=quantile(S2.prop.CIW, 0.25, na.rm=TRUE), 
              upper_S2=quantile(S2.prop.CIW, 0.75, na.rm=TRUE), median_S3=median(S3.prop.CIW, na.rm=TRUE), 
              lower_S3=quantile(S3.prop.CIW, 0.25, na.rm=TRUE), upper_S3=quantile(S3.prop.CIW, 0.75, na.rm=TRUE))
  Prop.CIW2 <- Prop.CIW %>% pivot_longer(cols=starts_with(c("median_","upper_","lower_")), 
                                         names_to=c("value_type","Method"), values_to="values", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="values")

  Prop.CIW.AJE <- full_join(AJE_base, Prop.CIW2 %>% dplyr::select(Time,Correlation,Method,median,upper,lower))
  Prop.CIW.AJE$Vals <- ifelse(is.na(Prop.CIW.AJE$median),NA,
                              paste0(format(round(Prop.CIW.AJE$median,3),nsmall=3, trim=FALSE)," (",
                                     format(round(Prop.CIW.AJE$lower,3),nsmall=3, trim=TRUE),"\U2013",
                                     format(round(Prop.CIW.AJE$upper,3),nsmall=3, trim=TRUE),")"))
  Prop.CIW.AJE$star <- ifelse(Prop.CIW.AJE$upper > 0.5,1,0)
  
  p.2.top <- ggplot(data=Prop.CIW.AJE, aes(x=median, y=Y)) +
    geom_point(size=2.5) + geom_errorbarh(mapping=aes(xmin=lower, xmax=ifelse(upper>0.5,NA,upper), y=Y), 
                                          inherit.aes=FALSE) +
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=8, color="black"),
          axis.line.x=element_blank(),
          strip.text=element_text(size=8),
          axis.ticks.x=element_line(color="black", size=1/ggplot2:::.pt), 
          text=element_text(size=8, color="black"),
          plot.tag=element_text(size=8, color="black"),
          plot.title=element_blank()) +
    labs(tag=paste0(ifelse(icc==corrs[3],FigLab[4],FigLab[1]),")")) +
    scale_x_continuous(name="Median", limit=c(-0.62,1.08), ## Limit is big enough to allow the text labels
                       breaks=seq(0,0.5,by=.1), labels=as.character(seq(0,0.5,by=.1))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,5.6), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=-0.62, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=-.33, y=5.2, label='underline("Day and Sampling")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) + 
    geom_text(aes(x=.56, y=Y, label=Vals, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=.85, y=5.2, label='underline("Median (IQR)")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) +
    annotate("segment", x=-.00295, xend=.50295, y=0, yend=0, color="black", size=1) ## Bottom axis that ends at break labels
  if (sum(Prop.CIW.AJE$star, na.rm=TRUE)>0) {
    StarVals <- Prop.CIW.AJE %>% filter(star==1)
    p.2.top <- p.2.top + geom_segment(data=StarVals,
                                  mapping=aes(x=median,xend=lower,y=Y,yend=Y), 
                                  arrow=arrow(angle=90, length=unit(0.03, "native"))) +
      geom_segment(data=StarVals,
                   mapping=aes(x=median,xend=pmin(0.5,upper),y=Y,yend=Y), 
                   arrow=arrow(length=unit(0.03, "native")))
  }
  
  ggsave(filename=paste0(plot_wd,"/Individual_Panels/Fig",ifelse(icc==corrs[3],FigNo[3],FigNo[2]),"_",
                         ifelse(icc==corrs[3],FigLab[4],FigLab[1]),".eps"),
         plot=p.2.top, width=3.5, height=2.5, units="in", device=cairo_ps)
  
  assign(x=paste0("p.",ifelse(icc==corrs[3],"S","2"),".",
                  ifelse(icc==corrs[3],FigLab[4],FigLab[1])), value=p.2.top)
  
  Prop.CICov <- ResFull_k %>% group_by(Time, Correlation) %>% 
    summarize(NSims=n(), CovIn_S1=mean(S1.prop.cvgIn, na.rm=TRUE), CovIn_S2=mean(S2.prop.cvgIn, na.rm=TRUE),
              CovIn_S3=mean(S3.prop.cvgIn, na.rm=TRUE), CovOut_S1=mean(S1.prop.cvgOut, na.rm=TRUE),
              CovOut_S2=mean(S2.prop.cvgOut, na.rm=TRUE), CovOut_S3=mean(S3.prop.cvgOut, na.rm=TRUE))
  Prop.CICov2 <- Prop.CICov %>% pivot_longer(cols=starts_with("Cov"), names_to=c("value_type","Method"), 
                                             values_to="Coverage", names_sep="_") %>%
    pivot_wider(names_from="value_type", values_from="Coverage")

  Prop.Cov.AJE <- full_join(AJE_base, Prop.CICov2 %>% dplyr::select(Time,Correlation,Method,CovIn))
  Prop.Cov.AJE$Cov <- Prop.Cov.AJE$CovIn*100
  Prop.Cov.AJE$Vals <- ifelse(is.na(Prop.Cov.AJE$Cov),NA,
                              paste0(format(round(Prop.Cov.AJE$Cov,1),nsmall=1)))
  
  p.2.bot <- ggplot(data=Prop.Cov.AJE, aes(x=Cov, y=Y)) +
    geom_point(size=2.5) + 
    theme_minimal() + 
    theme(panel.grid=element_blank(), axis.title.y=element_blank(), 
          axis.text.y=element_blank(), axis.text.x=element_text(size=8, color="black"),
          axis.title.x=element_text(size=8, color="black", hjust=.54),
          axis.line.x=element_blank(),
          axis.ticks.x=element_line(color="black", size=1/ggplot2:::.pt), 
          text=element_text(size=8, color="black"),
          strip.text=element_text(size=8),
          plot.tag=element_text(size=8, color="black"),
          plot.title=element_blank()) +
    labs(tag=paste0(ifelse(icc==corrs[3],FigLab[5],FigLab[2]),")")) +
    scale_x_continuous(name="Coverage, %", limit=c(78,110), ## Limit is big enough to allow the text labels
                       breaks=seq(90,100,by=5), labels=as.character(seq(90,100,by=5))) + ## Breaks & labels only where there is actual data
    scale_y_continuous(limit=c(0,5.6), expand=expansion(add=c(0,0))) + ## Needs to match Y vector, but with some space above and below
    geom_text(aes(x=78, y=Y, label=label, hjust="left"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Left side text labels
    annotate("text", x=83.5, y=5.2, label='underline("Day and Sampling")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) + 
    geom_text(aes(x=105, y=Y, label=Vals, hjust="center"), inherit.aes=FALSE, size=8/ggplot2:::.pt) + ## Right side text labels
    annotate("text", x=105, y=5.2, label='underline("Coverage, %")', parse=TRUE, hjust="center", size=8/ggplot2:::.pt) +
    annotate("segment", x=89.9415, xend=100.0585, y=0, yend=0, color="black", size=1) + ## Bottom axis that ends at break labels
    geom_segment(x=95, xend=95, y=0, yend=4.8, size=1/ggplot2:::.pt)
  
  ggsave(filename=paste0(plot_wd,"/Individual_Panels/Fig",ifelse(icc==corrs[3],FigNo[3],FigNo[2]),"_",
                         ifelse(icc==corrs[3],FigLab[5],FigLab[2]),".eps"),
         plot=p.2.bot, width=3.5, height=2.5, units="in", device=cairo_ps)
  
  assign(x=paste0("p.",ifelse(icc==corrs[3],"S","2"),".",
                  ifelse(icc==corrs[3],FigLab[5],FigLab[2])), value=p.2.bot)
  }
  
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[1],".eps"),
         plot=p.1.A+p.1.B+p.1.C+p.1.D+p.1.E+p.1.F+plot_layout(nrow=3,ncol=2,byrow=TRUE),
         width=7, height=7.5, units="in", device=cairo_ps)
  
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[2],".eps"),
         plot=p.2.A+p.2.B+p.2.C+p.2.D+plot_layout(nrow=2,ncol=2,byrow=TRUE),
         width=7, height=5, units="in", device=cairo_ps)
  
  ggsave(filename=paste0(plot_wd,"/Fig",FigNo[3],".eps"),
         plot=p.S.A+p.S.B+p.S.C+p.S.D+p.S.E+plot_layout(nrow=3,ncol=2,byrow=TRUE),
         width=7, height=7.5, units="in", device=cairo_ps)
}

## Proportion of Simulations with at least one symptomatic case, by analysis method
Prop.Obs <- ResFull %>% group_by(Time,Correlation,k) %>% 
  summarize(Prop.S1=mean(S1Symps>0), Prop.S2=mean(S2Symps>0), Prop.S3=mean(S3Symps>0))
Prop.Obs.t <- ResFull %>% group_by(Time) %>% 
  summarize(Prop.S1=mean(S1Symps>0), Prop.S2=mean(S2Symps>0), Prop.S3=mean(S3Symps>0))


