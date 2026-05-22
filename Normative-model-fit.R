library(readxl)   
library(gamlss)   
library(dplyr)    
library(doParallel)
library(foreach)    
library(iterators) 
library(reshape2) 
library(ggplot2)

datapath='D:/A_Normative_model/code_and_tables/code/'
setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

clinical_datapath=paste0(datapath,'20260405_update.csv')

savepath='D:/A_Normative_model/code_and_tables/code/Results/Step1_normative_models111';


MR_datapath=paste0(datapath,'MR_Results_limbic_20260404.xlsx')
var<-c('limb.vol.table')

library(stringr)

for(sheet in var) 
{ 
  
  
  setwd(datapath)
  MRI <- read_excel(MR_datapath,sheet=sheet)
  
  MRI<-MRI[!is.na(MRI$Freesurfer_Path3),]
  MRI<-as.data.frame(MRI)
  
  rownames(MRI)<-paste0(MRI$Freesurfer_Path2,MRI$Freesurfer_Path3)
  
  if(str_detect(sheet,'limb'))
  {
    tem_feature<-colnames(MRI)[c(2:15)];
  }
    
    str=sheet;
  
  
  if (!(dir.exists(paste0(savepath,'/',str))))
  {dir.create(paste0(savepath,'/',str))}
  
  setwd(paste0(savepath,'/',str))
  
  
  setwd(datapath)
  data1<-read.csv(clinical_datapath,header=TRUE);
  data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)
  
  
  Quant_data<-list()
  
  
  for(i in tem_feature[1:length(tem_feature)])
  {
    
    print(i)
    setwd(paste0(savepath,'/',str))
    if(file.exists(paste0(str,'_',i,'_loop_our_model.rds'))){print('file exist');next;print('file exist')};
    
    
    # for each feature, we should load the original clinical information
    setwd(datapath)  
    data1<-read.csv(clinical_datapath,header=TRUE);
    data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)
    
    site_count <- data1 %>%
      group_by(Site_ZZZ) %>%  
      summarise(count = n())  
    site_count<-site_count[order(site_count$count),]
    print(site_count)
    
    for(site in unique(site_count$Site_ZZZ))
    {
      if(site_count[site_count$Site_ZZZ==site,'count']<10)
      {
        data1[data1$Site_ZZZ==site,'Database_included']<-0
      }
      
    }
    
    rownames(data1)<-paste0(data1$Freesufer_Path2,data1$Freesufer_Path3)
    
    setwd(paste0(savepath,'/',str))
    
    
    inter_row<-intersect(rownames(data1),rownames(MRI))
    data1=cbind(data1[inter_row,],MRI[inter_row,i])
    
    
    colnames(data1)[dim(data1)[2]]=c('tem_feature')
    
    all_data<-data1[data1$Database_included==1&
                      !is.na(data1$baseline),]
    
    rownames(all_data)<-paste0(all_data$Freesufer_Path2,all_data$Freesufer_Path3)
    all_data_original<-all_data
    
    data1<-all_data
    
    data1=data1[!is.na(data1$tem_feature)&!is.na(data1$Data_baseline)&data1$Diagnosis=='HC'&data1$Age>=4&data1$Age<=85,]
    
    data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    
    data1$Sex<-as.factor(data1$Sex)
    
    data1$Sex<-factor(data1$Sex,levels=c('Female','Male'))
    
    
    data1<-data1[order(data1$Age),]
    data1[,'feature']<-data1$tem_feature
    all_data[,'feature']<-all_data$tem_feature
    
    #remove the extreme values
    data1<-data1[!is.na(data1$tem_feature),]
    data1<-data1[data1$feature>(mean(data1$feature)-3*sd(data1$feature))&
                   data1$feature<(mean(data1$feature)+3*sd(data1$feature)),]
    
    data1<-data1[data1$feature>0,]
    #select the columns
    data1<-data1[,c('Age','Sex','Site_ZZZ','tem_feature','feature')]
    
    data1_backup<-data1;
    data1_child<-data1[data1$Age<=18,];
    data1_adult<-data1[data1$Age>18&data1$Age<70,];
    data1_old<-data1[data1$Age>=70,];
    data1_adult_sample<- data1_adult %>% sample_frac(0.3)
    data1<-rbind(data1_child,data1_adult_sample,data1_old)
    
    
    list_par<-data.frame(matrix(0,3*3*2*2,1));
    
    con=gamlss.control()
    
    num=0;
    
    results_try<-try({
      for(i_poly in 1:3)
      {
        for(j_poly in 1:3)
        {
          for(i_rnd in 0:1)
          {
            for(j_rnd in 0:1)
            {
              num=num+1;
              list_par[num,1]<-i_poly
              list_par[num,2]<-j_poly
              list_par[num,3]<-i_rnd
              list_par[num,4]<-j_rnd
              
            }}}}
      library(doParallel)
      library(foreach)
      
      
      cl<-makeCluster(10)
      registerDoParallel(cl)
      my_data<-foreach(num=1:dim(list_par)[1],
                       .combine=rbind,
                       .packages = c('gamlss')) %dopar% fit_model(num)
      stopCluster(cl)
      
      
      list_fit<-my_data
      print(list_fit)
      
      
      #fit using the bestfit npoly and random with lowest BIC
      model_ind<-which.min(list_fit$BIC);
      sel_mu_poly=list_fit$mu_poly[model_ind]
      sel_sigma_poly=list_fit$sigma_poly[model_ind]
      i_rnd=list_fit$mu_random[model_ind]
      j_rnd=list_fit$sigma_random[model_ind]
    })
    
    
    
    if(inherits(results_try,'try-error')) 
    {sel_mu_poly=2
    sel_sigma_poly=2
    i_rnd=1
    j_rnd=1
    con=gamlss.control(c.crit = 0.01, n.cyc = 2,autostep = FALSE)
    }
    
    
    
    data1<-data1_backup;
    
    
    m0<-best_fit(sel_mu_poly,sel_sigma_poly,i_rnd,j_rnd)
    
    
    if(i_rnd==1&j_rnd==1){
      m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data1)}else if(i_rnd==1&j_rnd==0){
                   m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                              control=con,
                              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                              data=data1)}else if(i_rnd==0&j_rnd==1){
                                m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                                           control=con,
                                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                           data=data1)}else if(i_rnd==0&j_rnd==0){
                                             m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                                                        control=con,
                                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                        data=data1)}
    
    
    #for all population plot
    if(i_rnd==1&j_rnd==1){
      m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
                 sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data1)}else if(i_rnd==1&j_rnd==0){
                   m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
                              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
                              control=con,
                              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                              data=data1)}else if(i_rnd==0&j_rnd==1){
                                m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
                                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
                                           control=con,
                                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                           data=data1)}else if(i_rnd==0&j_rnd==0){
                                             m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
                                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
                                                        control=con,
                                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                        data=data1)}
    model1<-m3;
    
    num_length=5000
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
    tem_par<-mu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    mu=par
    
    
    tem_par<-sigma0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    sigma=par
    
    
    tem_par<-nu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    nu=par
    
    
    p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                 cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                 calibration=FALSE,lpar=3)
    p2[,'sigma']<-sigma
    
    
    
    
    
    model1<-m2;
    
    
    
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
    tem_par<-mu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    mu=par
    
    
    tem_par<-sigma0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    sigma=par
    
    tem_par<-nu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    nu=par
    
    
    p2_all<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                     cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                     calibration=FALSE,lpar=3)
    p2_all[,'sigma']<-sigma
    
    
    
    #male
    model1<-m2;
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
    tem_par<-mu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    mu=par
    
    
    tem_par<-sigma0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    sigma=par
    
    tem_par<-nu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    nu=par
    
    
    male_p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                      cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                      calibration=FALSE,lpar=3)
    male_p2[,'sigma']<-sigma
    
    
    #female
    
    model1<-m2;
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
    tem_par<-mu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    mu=par
    
    
    tem_par<-sigma0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    sigma=par
    
    tem_par<-nu0
    par<-tem_par[1:num_length]
    Seg=length(tem_par)/num_length;
    if(Seg>1)
    {
      for(Seg1 in c(2:Seg))
      {
        par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
      }
      par=par/Seg
    }
    nu=par
    
    female_p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                        cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                        calibration=FALSE,lpar=3)
    female_p2[,'sigma']<-sigma
    
    
    
    library(reshape2);
    colnames(p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    mydata<-melt(p2,id='Age');colnames(mydata)<-c('Age','Percentile','Value')
    
    step_age<-(max(data1$Age)-min(data1$Age))/num_length
    dim(p2)[1]-1
    Grad_p2<-(p2$median[2:dim(p2)[1]]-p2$median[1:(dim(p2)[1]-1)])/step_age
    Grad_p2<-data.frame(c(Grad_p2,Grad_p2[dim(p2)[1]-1]));
    p2<-cbind(p2,Grad_p2)
    colnames(p2)[dim(p2)[2]]<-c('Gradient1')
    
    if(!(str_detect(sheet,'thickness')))
    {
      scale1=10000;
      ylab1='×10^4 mm3';
    }
    
    if(str_detect(sheet,'thickness'))
    {
      scale1=1;
      ylab1='mm';
    }
    
    png(filename = paste0(str,'_',i,'_all_without_sex_stratified_Gradient.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      geom_line(data=p2,aes(x=Age,y=Gradient1/scale1),color=c('#262626'),linewidth=1,linetype=c('solid'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",        
          #face = "bold" ,         
          family = "serif"
        )
      )+
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    
    print(p3)  
    dev.off()
    
    
    
    png(filename = paste0(str,'_',i,'_all_without_sex_stratified.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#E84935'),shape=16,size=3,alpha = 0.1)+
      geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#4FBBD8'),shape=17,size=3,alpha = 0.1)+
      geom_line(data=p2,aes(x=Age,y=median/scale1),color=c('#262626'),linewidth=1,linetype=c('solid'))+
      geom_line(data=p2,aes(x=Age,y=lower99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
      geom_line(data=p2,aes(x=Age,y=lower95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=p2,aes(x=Age,y=upper95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=p2,aes(x=Age,y=upper99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"     
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",         
          family = "serif"
        )
      )+
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    dev.off()
    
    
    #for all population with sex stratified
    #female data
    colnames(female_p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    mydata<-melt(female_p2,id='Age');colnames(mydata)<-c('Age','Percentile','Value')
    Female_p2<-female_p2;
    
    
    #male data
    colnames(male_p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    mydata<-melt(male_p2,id='Age');colnames(mydata)<-c('Age','Percentile','Value')
    Male_p2<-male_p2;
    
    
    dim(Female_p2)[1]-1
    Grad_Female_p2<-(Female_p2$median[2:dim(Female_p2)[1]]-Female_p2$median[1:(dim(Female_p2)[1]-1)])/step_age
    Grad_Female_p2<-data.frame(c(Grad_Female_p2,Grad_Female_p2[dim(Female_p2)[1]-1]));
    Female_p2<-cbind(Female_p2,Grad_Female_p2)
    colnames(Female_p2)[dim(Female_p2)[2]]<-c('Gradient1')
    
    dim(Male_p2)[1]-1
    Grad_Male_p2<-(Male_p2$median[2:dim(Male_p2)[1]]-Male_p2$median[1:(dim(Male_p2)[1]-1)])/step_age
    Grad_Male_p2<-data.frame(c(Grad_Male_p2,Grad_Male_p2[dim(Male_p2)[1]-1]));
    Male_p2<-cbind(Male_p2,Grad_Male_p2)
    colnames(Male_p2)[dim(Male_p2)[2]]<-c('Gradient1')
    
    
    png(filename = paste0(str,'_',i,'_all_with_sex_stratified_Gradient.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      geom_line(data=Female_p2,aes(x=Age,y=Gradient1/scale1),color=c('#E84935'),linewidth=1,linetype=c('solid'))+
      geom_line(data=Male_p2,aes(x=Age,y=Gradient1/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('solid'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",         
          family = "serif"
        )
      )+
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    dev.off()
    
    
    png(filename = paste0(str,'_',i,'_all_with_sex_stratified_sigma.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      geom_line(data=Female_p2,aes(x=Age,y=sigma),color=c('#E84935'),linewidth=1,linetype=c('solid'))+
      geom_line(data=Male_p2,aes(x=Age,y=sigma),color=c('#4FBBD8'),linewidth=1,linetype=c('solid'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",         
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",         
          family = "serif"
        )
      )+
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    dev.off()
    
    
    
    
    png(filename = paste0(str,'_',i,'_all_with_sex_stratified.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#E84935'),shape=16,size=3,alpha = 0.1)+
      geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#4FBBD8'),shape=17,size=3,alpha = 0.1)+
      geom_line(data=Female_p2,aes(x=Age,y=median/scale1),color=c('#E84935'),linewidth=1,linetype=c('solid'))+
      geom_line(data=Female_p2,aes(x=Age,y=lower95CI/scale1),color=c('#E84935'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=Female_p2,aes(x=Age,y=upper95CI/scale1),color=c('#E84935'),linewidth=1,linetype=c('dotted'))+
      
      geom_line(data=Male_p2,aes(x=Age,y=median/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('solid'))+
      geom_line(data=Male_p2,aes(x=Age,y=lower95CI/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=Male_p2,aes(x=Age,y=upper95CI/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('dotted'))+
      
      
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,             
          color = "black",         
          family = "serif"
        )
      )+
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    
    dev.off()
    
    #calculate the quantile for all cases including HC and Diseases
    
    Quant_score_sum<-NULL
    
    
    all_data<-all_data[all_data$tem_feature!=''&!is.null(all_data$tem_feature)&!is.na(all_data$tem_feature)&!is.infinite(all_data$tem_feature),]
    all_data<-all_data[all_data$Age!=''&!is.null(all_data$Age)&!is.na(all_data$Age)&!is.infinite(all_data$Age),]
    all_data<-all_data[all_data$Site_ZZZ!=''&!is.null(all_data$Site_ZZZ)&!is.na(all_data$Site_ZZZ)&!is.infinite(all_data$Site_ZZZ),]
    all_data<-all_data[all_data$Sex!=''&!is.null(all_data$Sex)&!is.na(all_data$Sex)&!is.infinite(all_data$Sex),]
    
    all_data1<-all_data[,c('Age','Sex','Site_ZZZ','tem_feature')]
    
    all_data1$Sex<-as.factor(all_data1$Sex)
    all_data1$Site_ZZZ<-as.factor(all_data1$Site_ZZZ)
    
    
    model1<-m2;
    
    
    for(sub in 1:dim(all_data)[1])
    {
      
      if(!is.null(m2$mu.coefSmo[[1]]))
      {
        if(!(all_data1$Site_ZZZ[sub] %in% names(m2$mu.coefSmo[[1]]$coef)))
        {
          all_data1$Site_ZZZ[sub]<-names(which.max(abs(m2$mu.coefSmo[[1]]$coef-mean(m2$mu.coefSmo[[1]]$coef))))
          
        }
      }
      
      if(!is.null(m2$sigma.coefSmo[[1]]))
      {
        if(!(all_data1$Site_ZZZ[sub] %in% names(m2$sigma.coefSmo[[1]]$coef)))
        {
          all_data1$Site_ZZZ[sub]<-names(which.max(abs(m2$sigma.coefSmo[[1]]$coef-mean(m2$sigma.coefSmo[[1]]$coef))))
        }
      }
      
    }
    
    
    mu <- predict(model1, newdata = all_data1, type = "response", what = "mu")
    sigma <- predict(model1, newdata = all_data1, type = "response", what = "sigma")
    nu <- predict(model1, newdata = all_data1, type = "response", what = "nu")
    
    
    
    if(length(mu)!=dim(all_data1)[1])
    {
      print("Error, Please Check Data!")
    }
    
    
    Quant_score_sum<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
                              xname = 'Age',xvalues=all_data1$Age,yval=all_data1$tem_feature,
                              calibration=FALSE,lpar=3,cdf=TRUE)
    
    
    Quant_score_sum<-data.frame(Quant_score_sum);
    colnames(Quant_score_sum)<-c('Quant_score');
    rownames(Quant_score_sum)<-rownames(all_data1)
    
    Quant_data[[i]]<-Quant_score_sum
    
    
    results<-list();
    results$Female_p2<-Female_p2
    results$Male_p2<-Male_p2
    results$p2<-p2
    results$peakage<-p2$Age[which.max(p2$median)]
    results$p2_all<-p2_all
    results$m2<-m2
    results$m0<-m0
    results$m3<-m3
    
    results$list_fit<-list_fit
    
    results$Quant_data<-Quant_data
    results$data1<-data1
    results$all_data<-all_data1
    results$str<-str
    results$i<-i
    results$all_data_original<-all_data_original
    
    
    saveRDS(results,paste0(str,'_',i,'_loop_our_model.rds'))
    
  }
  
}



