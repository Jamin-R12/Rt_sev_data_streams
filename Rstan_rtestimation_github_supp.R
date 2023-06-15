
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
#3) RT estimation via cmdstan ####
dv_case<-AggCases %>%
  mutate(mean_rt = NA)%>%
  arrange(date_int) %>%
  mutate(day = row_number())

#Set variables for CMDStan model
dv_case$newI_case<-round(as.numeric(dv_case$newI_case))
  
T<-nrow(dv_case) #T, number of event time
S<-12 #S, Serial interval 
I<-dv_case$newI_case %>%as.vector()    
w<-c(0.000,0.036, 0.167, 0.208, 0.180, 0.135, 0.096, 0.066, 0.045, 0.031, 0.021, 0.015) 
tau<- 7 
Est_rt_stan<-list(T=T, S=S, I=I, w=w, tau=tau)

model<-cmdstan_model('temp_data/rt_est_7_ex.stan')

fit<-model$sample(
  data = Est_rt_stan,
  seed = 123456,
  chains = 4,
  parallel_chains = 4,
  refresh = 500)
  
stanfit <-read_stan_csv(fit$output_files())

rt_summary <- summary(stanfit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything())
  
#change data. 
rt_case<-rt_summary %>% 
  mutate(
    mean_rt = mean,
    lower = rt_summary$`2.5%`,
    upper = rt_summary$`97.5%`,
    day = row_number())%>%
  select(day, mean_rt, lower, upper, sd, se_mean, n_eff)%>%
  filter(!abs(mean_rt)>(100*mean(mean_rt)))


#plot
plot(rt_case$day, rt_case$mean_rt, type = "l", xlab = "Day", ylab = "Rt")
lines(rt_case$day, rt_case$upper, type = "l", lwd = 2, col = "grey")
lines(rt_case$day, rt_case$lower, type = "l", lwd = 2, col = "grey")
lines(rt_case$day, rt_case$mean_rt, type = "l", lwd = 2)
legend("topright", legend = c("Instantaneous Rt", "Credible interval"), 
       lty = c(1, 1), lwd = c(2, 2), col = c("black", "grey"))













