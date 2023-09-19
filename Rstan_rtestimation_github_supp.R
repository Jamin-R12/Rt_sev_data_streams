
#Introduction#####
#script shows and example using aggregated case counts, it includes the joint PMF calculation, deconvolution and rt estimation.  

#packages
pacman::p_load(dplyr, readr,  rstan, cmdstanr, posterior, bayesplot, surveillance)


#----------------------------------------------------------------------#
#1) Joint PMF####
OnsetToConfirmed<-read_csv("output/OnsetToConfirmed.csv")

#onset to confirmed PMF
time_conf<-OnsetToConfirmed%>%
  filter(!is.na(conf_ons) & conf_ons>=-2 & conf_ons<=14)%>%
  mutate(N = n())%>%
  group_by(conf_ons)%>%
  summarise(dist = n(),
            pmf  = n()/N)%>%
  slice_head(n = 1)
conf<-list(day=-2:14, prob=time_conf$pmf)

#incubation period PMF
mean = 3.1
sd = 2.3
m.log<-log(mean)-0.5*(log((sd/mean)^2+1)) #from where???
sd.log<-sqrt(log((sd/mean)^2 +1))
incub<-list(day=0:14,prob=dlnorm(0:14,m.log,sd.log))

j.incub.conf<-matrix(nrow=(15+length(conf$day)),ncol = 2)
for (i in 1:length(incub$day)){
  for (j in 1:length(conf$day)){
    if (is.na(j.incub.conf[i+j,2])){                    # if first prob. if missing already
      
      j.incub.conf[i+j,1]=(incub$day[i]+conf$day[j])
      j.incub.conf[i+j,2]=(incub$prob[i]*conf$prob[j])
      
    }
    else if (!is.na(j.incub.conf[i+j,2])){                       # if following prob. add in
      
      j.incub.conf[i+j,1]=(incub$day[i]+conf$day[j])
      j.incub.conf[i+j,2]=j.incub.conf[i+j,2]+(incub$prob[i]*conf$prob[j])
    }
    if(j.incub.conf[i+j,1]<=0){
      j.incub.conf[i+j,2]<-0
    } 
  }
}
j.incub.conf<-j.incub.conf[-(1:3),]

plot(j.incub.conf[,1],j.incub.conf[,2])

j.pmf.conf<-j.incub.conf[,2]/(sum(j.incub.conf[,2]))





#----------------------------------------------------------------------#
#2) Deconvolution####

dcv<-function(dv_set, indi){
  print(dv_set)
  ZZ<-sts(dv_set[,paste0(as.character(indi),"_cnt")])
  bpnp.control <- list(k=0,
                       eps=rep(1,2),
                       iter.max=rep(250,2),
                       B=1,               
                       verbose=TRUE)     
  set.seed(123)
  temp <- backprojNP(ZZ,
                     incu.pmf=j.pmf.conf, 
                     control=modifyList(bpnp.control,
                                        list(eq3a.method = "C")
                     ),
                     ylim=c(0,max(X,Y)))
  newI2 <<- temp@upperbound
}

##read-in case aggregate data. 
AggCases<-read_csv(paste0("output/Aggregate_case_wave5") )

#run back projectin. 
dcv(AggCases,"case")
newI2 <- data.frame(matrix(unlist(newI2), nrow=length(newI2), byrow=TRUE))

AggCases$newI_case<-newI2[,1]
AggCases[is.na(AggCases)] <- 0
AggCases<-AggCases%>%
  filter(!newI_case==0)

#plot deconvoluted datastream. 
AggCases %>%
  ggplot(aes(x=date_int))+
  geom_line(aes(y=newI_case),colour= "red" )+
  geom_line(aes(y=case_cnt), colour = "blue")




#----------------------------------------------------------------------#
#3) RT estimation via EpiEstim ####
dv_case<-AggCases 

  dv_case<-dv_case %>% 
    filter(date_int>=START_DATE)
  
  T_row <- nrow(dv_case)
  t_start <- seq(2, T_row-6)
  t_end <- t_start + 6 
  res_para_temp <- estimate_R(incid = dv_case$newI, 
                              method="non_parametric_si",
                              config = make_config(list(
                                t_start = t_start,
                                t_end = t_end,
                                si_distr = c(0.000,0.036, 0.167, 0.208, 0.180, 0.135, 0.096, 0.066, 0.045, 0.031, 0.021, 0.015))
                              ))
  

  
  rt_case = list(Estim = res_para_temp,
               df = as.data.frame(cbind(date = dv_case[t_end,'date_int'],
                                        res_para_temp$R,
                                        indi = 'case',
                                        wave = WAVE))
  )

res_L
#plot
plot(rt_case$df$date_int,rt_case$df$`Mean(R)`, type = "l", xlab = "Day", ylab = "Rt")
lines(rt_case$df$date_int, icu_L_6$df$`Quantile.0.025(R)`, type = "l", lwd = 2, col = "grey")
lines(rt_case$df$date_int, icu_L_6$df$`Quantile.0.975(R)`, type = "l", lwd = 2, col = "grey")
legend("topright", legend = c("Instantaneous Rt", "Credible interval"), 
       lty = c(1, 1), lwd = c(2, 2), col = c("black", "grey"))













