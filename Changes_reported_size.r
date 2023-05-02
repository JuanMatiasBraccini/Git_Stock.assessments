
  Logbook=Logbook%>%
    filter(LatDeg<=(-26) & method=='GN')%>%
    mutate(species=ifelse(species==19000 & LongDeg>116,19004,species))%>%
    group_by(Same.return.SNo,species,finyear,mshigh)%>%
    summarise(nfish=sum(nfish),
              livewt=sum(livewt))%>%
    ungroup()%>%
    dplyr::select(finyear,livewt,nfish,species,mshigh)%>%
    mutate(Mean.wght=livewt/nfish)%>%
    left_join(All.species.names,by=c("species"="SPECIES"))%>%
    filter(SNAME%in%names(Species.data))%>%
    left_join(Wei.range%>%dplyr::select(TW.min,TW.max,SPECIES),by=c("species"="SPECIES"))%>%
    filter(Mean.wght>=TW.min & Mean.wght<=TW.max)
  
  N.min=Logbook%>%
    group_by(finyear,species,mshigh)%>%
    tally()%>%
    filter(n>=Min.annual.obs.ktch)%>%
    mutate(Keep="YES")
  
  Logbook=Logbook%>%
    left_join(N.min,by=c('finyear','species','mshigh'))%>%
    filter(Keep=="YES")%>%
    mutate(Finyear=as.numeric(substr(finyear,1,4)),
           Finyear.d=factor(Finyear,levels=sort(unique(Finyear))))%>%
    filter(!is.na(Mean.wght))%>%
    mutate(Mesh=ifelse(mshigh==165,"6.5",ifelse(mshigh==178,'7',NA)))%>%
    filter(!is.na(Mesh))
  
  Logbook.sp=sort(unique(Logbook$SNAME))
  Logbook.sp=subset(Logbook.sp,!Logbook.sp=="angel sharks")
  
  Logbook=Logbook%>%filter(Finyear<=as.numeric(substr(Last.yr.ktch,1,4))) #keep only years with catch 
  
  Change.mean.weight.catch=vector('list',length(Logbook.sp))
  names(Change.mean.weight.catch)=Logbook.sp
  any_column_NA <- function(x) any(is.numeric(x))
  replace_pos <- function(x) if_else(x>0,1,x)
  fn.plt.mn.ktch.wght=function(d.list,NM,XLIM,show.data=FALSE,wei.mat,max.wei,recruit.cpue)
  {
    my_formula = y ~ x
    if(nrow(d.list)>0 & length(unique(d.list$finyear))>=2)
    {
      Kip=d.list%>%
        group_by(finyear,mshigh)%>%
        tally()%>%
        spread(mshigh,n,fill=0)%>%
        mutate_if(any_column_NA,replace_pos)%>%
        gather(mesh,n,-finyear)%>%
        group_by(mesh)%>%
        summarise(n=sum(n))%>%
        filter(n>1)
      
      p=d.list%>%
        filter(mshigh%in%as.numeric(unique(Kip$mesh)))%>%
        mutate(Mesh=paste(Mesh,"inch"))%>%
        ggplot(aes(x = Finyear, y = Mean.wght)) 
      if(show.data) p=p+geom_point(position = position_jitter(seed = 1, width = 0.2),color='grey',alpha=.2)
      p=p+
        geom_violin(aes(fill = Finyear, group = Finyear.d,color=Finyear), alpha = 0.2) + 
        facet_wrap(~Mesh,scales='free')+
        stat_summary(fun = "mean",geom = "point",color = "red",size=2)+
        geom_smooth(color = rgb(1,.1,.1,alpha=0.2), formula = my_formula, method = 'lm',se=TRUE)+ 
        stat_poly_eq(aes(label = paste("atop(", stat(eq.label),  ",", 
                                       paste(stat(adj.rr.label),stat(p.value.label), sep = "*\", \"*"), ")")),
                     formula = my_formula, parse = TRUE,
                     label.y = "top", label.x = "right", size = 4.5) +
        xlab("")+ylab("")+
        ggtitle(NM)+xlim(XLIM)+
        theme_PA(axs.T.siz=22,axs.t.siz=13,str.siz=16)+
        theme(legend.position = "none",
              plot.title =element_text(size=17))+
        ylim(0,max.wei)+
        geom_hline(yintercept = max.wei,colour='orange',size=1.05)+
        geom_hline(yintercept = wei.mat,colour='orange',linetype = "dashed",size=1.05)
      
      return(p)
    }
    
  }
  for(s in 1:length(Logbook.sp)) 
  {
    print(paste("Change in mean weight of landed individual ","--",Logbook.sp[s]))
    id=match(Logbook.sp[s],names(List.sp))
    dummy=print(fn.plt.mn.ktch.wght(d.list=Logbook%>%filter(SNAME==Logbook.sp[s]),
                                    NM=capitalize(Logbook.sp[s]),
                                    XLIM=c(min(Logbook$Finyear)-1,1+max(Logbook$Finyear)),
                                    wei.mat=with(List.sp[[id]],AwT*TL.50.mat^BwT),  
                                    max.wei=with(List.sp[[id]],AwT*TLmax^BwT),
                                    recruit.cpue=recruit.cpiui[[id]]))
    if(!is.null(dummy)) Change.mean.weight.catch[[s]]=dummy
    rm(dummy)
  }
  
  #Get catch rates of recruits
  fun.get.prop.zero.plus.length=function(d,NM,toMatch,min.annual.obs,size.1.yr.old)
  {
    d.list=d[grep(paste(toMatch,collapse="|"),names(d))]
    if(sum(grepl('Table',names(d.list)))>0)d.list=d.list[-grep('Table',names(d.list))]
    
    if(length(d.list)>0)
    {
      #by myear
      d.list=do.call(rbind,d.list)
      d.list=d.list%>%
        mutate(mesh=case_when(grepl("6.5.inch",rownames(d.list))~6.5,
                              grepl("7.inch",rownames(d.list))~7),
               Zeroo.plus=ifelse(FL<=size.1.yr.old,'Yes','No'))
      N.min=d.list%>%
        group_by(FINYEAR,Zeroo.plus)%>%  
        tally()%>%
        spread(Zeroo.plus,n,fill=0)%>%
        mutate(n=No+Yes)%>%
        filter(n>=min.annual.obs)%>%
        mutate(Prop.zero.plus=Yes/n)
      if(nrow(N.min)>0)return(N.min%>%dplyr::select(FINYEAR,n,Prop.zero.plus))
    }
  }
  fun.get.prop.zero.plus.wght=function(d.list,NM,min.annual.obs,whgt.1.yr.old)
  {
    if(nrow(d.list)>0 & length(unique(d.list$finyear))>=2)
    {
      Kip=d.list%>%
        filter(nfish==1)%>%    #if more than 1 reported shark, it could be one big and one small for e.g. so don't know actual individual weight     
        mutate(Zeroo.plus=ifelse(Mean.wght<=whgt.1.yr.old,'Yes','No'))
      
      N.min=Kip%>%
        group_by(finyear,Zeroo.plus)%>%  
        tally()%>%
        spread(Zeroo.plus,n,fill=0)
      if(!'Yes'%in%colnames(N.min)) N.min$Yes=0
      N.min=N.min%>%
        mutate(n=No+Yes,
               Prop.zero.plus=Yes/n)%>%
        rename(FINYEAR=finyear)
      # N.min=N.min%>%filter(n>=min.annual.obs)
      if(nrow(N.min)>0)return(N.min%>%dplyr::select(FINYEAR,n,Prop.zero.plus))
    }
  }
  recruit.cpiui=vector('list',N.sp)     
  names(recruit.cpiui)=Keep.species
  for(i in 1:length(recruit.cpiui))
  {
    #Get TDGDLF cpues
    CPUE=compact(Catch.rate.series[[i]])
    Kip=grep('TDGDLF',names(CPUE)) 
    CPUE=CPUE[Kip]
    DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
    if(length(DROP)>0)CPUE=CPUE[-DROP]
    
    if(length(CPUE)>0)
    {
      #Get proportion of 0+ individuals
      #Observed size frequency (not used, too short a series)
      prop.zero.length=fun.get.prop.zero.plus.length(d=Species.data[[i]],
                                                     NM=capitalize(names(Species.data)[i]),
                                                     toMatch=c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2"),
                                                     min.annual.obs=Min.annual.obs,
                                                     size.1.yr.old=ceiling(with(List.sp[[i]],Lzero+(Growth.F$FL_inf-Lzero)*(1-exp(-Growth.F$k*1)))))
      
      if(!is.null(prop.zero.length))
      {
        CPUE.prop.zero.plus.length=CPUE
        for(x in 1:length(CPUE))
        {
          CPUE.prop.zero.plus.length[[x]]=CPUE.prop.zero.plus.length[[x]]%>%
            left_join(prop.zero.length,by=c('Finyear'='FINYEAR'))%>%
            mutate(Mean.prop=Mean*Prop.zero.plus,
                   LOW.CI.prop=LOW.CI*Prop.zero.plus,
                   UP.CI.prop=UP.CI*Prop.zero.plus)
        }
      }
      
      #Reported weight (only applicable to daily)
      id=match(names(List.sp)[i],Logbook.sp)
      prop.zero.wght=fun.get.prop.zero.plus.wght(d.list=Logbook%>%filter(SNAME==Logbook.sp[id]),
                                                 NM=capitalize(Logbook.sp[id]),
                                                 min.annual.obs=Min.annual.obs,
                                                 whgt.1.yr.old=with(List.sp[[i]],AwT*(ceiling(Lzero+(Growth.F$FL_inf-Lzero)*(1-exp(-Growth.F$k*1)))*a_FL.to.TL+b_FL.to.TL)^BwT))
      CPUE.prop.zero.plus.wght=CPUE[["TDGDLF.daily"]]
      if(!is.null(prop.zero.wght) & !is.null(CPUE.prop.zero.plus.wght))
      {
        CPUE.prop.zero.plus.wght=CPUE.prop.zero.plus.wght%>%
          left_join(prop.zero.wght,by=c('Finyear'='FINYEAR'))%>%
          mutate(Mean.prop=Mean*Prop.zero.plus,
                 LOW.CI.prop=LOW.CI*Prop.zero.plus,
                 UP.CI.prop=UP.CI*Prop.zero.plus)
        Full.d=CPUE.prop.zero.plus.wght%>%
          dplyr::select(yr.f,Mean,LOW.CI,UP.CI)%>%
          mutate(dat='All ages')
        Rec.d=CPUE.prop.zero.plus.wght%>%
          dplyr::select(yr.f,Mean.prop)%>%
          mutate(LOW.CI=NA,
                 UP.CI=NA,
                 dat='0+')%>%
          rename(Mean=Mean.prop)
        
        p=rbind(Full.d,Rec.d)%>%
          mutate(Mesh='6.5 & 7 inch combined')%>%
          ggplot(aes(yr.f,Mean,color=dat))+
          geom_point(size=2)+
          geom_errorbar(aes(ymin=LOW.CI, ymax=UP.CI), width=.2)+
          ylim(0,max(CPUE.prop.zero.plus.wght$UP.CI))+
          theme_PA(axs.T.siz=22,axs.t.siz=14,str.siz=16)+
          ylab("")+xlab("")+ggtitle(capitalize(names(List.sp)[i]))+
          theme(legend.title=element_blank())+
          facet_wrap(~Mesh,scales='free')
        
        recruit.cpiui[[i]]=p
      }
    }
  }
  for(l in 1:length(Store.spatial.temporal.ktch))
  {
    figure2=ggarrange(plotlist=fun.find.in.list(x=recruit.cpiui[Lista.sp.outputs[[l]]]),
                      ncol=1,common.legend = TRUE)+
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
    FIG=ggarrange(plotlist=list(figure2),ncol=1)+
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
    annotate_figure(FIG,
                    left = text_grob("Relative catch rate", rot = 90,size=20,vjust=1),
                    bottom = text_grob("Financial year",size=20,vjust=-0.5))
    
    ggsave(paste(Rar.path,'/Relative catch rate_recruits_TDGDLF_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = 12.5,height = 10,compression = "lzw")
  }
  
  
  #combine mean weight and observed size
  for(l in 1:length(Store.spatial.temporal.ktch))
  {
    figure=ggarrange(plotlist=fun.find.in.list(x=Change.mean.weight.catch[Lista.sp.outputs[[l]]]),ncol=1)+
      theme(plot.margin = unit(c(-.5,0,0,-.5), "lines"))
    figure2=ggarrange(plotlist=fun.find.in.list(x=Change.mean.length[Lista.sp.outputs[[l]]]),ncol=1)+
      theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
    figure2=annotate_figure(figure2,
                            left= text_grob('Observed fork length (cm)', rot = 90,size=20,vjust=2))
    FIG=ggarrange(plotlist=list(figure,figure2),ncol=2)+
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
    annotate_figure(FIG,
                    left = text_grob("Total weight of caught individuals (kg)", rot = 90,size=20,vjust=1.25),
                    bottom = text_grob("Financial year",size=20,vjust=-1))
    ggsave(paste(Rar.path,'/Changes in size_TDGDLF_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = 12,height = 13,compression = "lzw")
  }
  
  do.model.based.mn.weight.ktch=FALSE
  if(do.model.based.mn.weight.ktch)
  {
    Mn.weit.ktch=lapply(Species.data, function(x) x[["annual.mean.size_relative"]])
    Mn.weit.ktch=fun.find.in.list(x=Mn.weit.ktch, Drop=names(Indicator.species))
    
    Change.mean.weight.catch=vector('list',length(Mn.weit.ktch))
    names(Change.mean.weight.catch)=names(Mn.weit.ktch)
    
    fn.plt.mn.ktch.wght=function(d,NM)
    {
      my_formula = y ~ x
      p=d%>%
        mutate(Finyear.f=factor(Finyear),
               Finyear=as.numeric(substr(Finyear,1,4)))%>%      
        ggplot(aes(x = Finyear, y = mean)) +
        geom_point(aes(colour=Finyear.f),size=3)+
        geom_smooth(color = "black", formula = my_formula, method = 'lm',se=F)+ 
        stat_poly_eq(aes(label =  paste(stat(eq.label),stat(adj.rr.label),stat(p.value.label),
                                        sep = "*\", \"*")),
                     formula = my_formula, parse = TRUE,
                     label.y = "bottom", label.x = "right", size = 5) +
        xlab("")+ylab("")+
        theme_PA(axs.T.siz=18,axs.t.siz=14,str.siz=16)+
        theme(legend.position = "none",
              plot.title =element_text(size=20),
              panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = "white"),
              panel.border = element_rect(colour = "black", fill=NA, size=1.25))+
        labs(title=NM)+ylim(0,max(d$mean,na.rm=T)) 
      return(p)
    }
    for(s in 1:length(Mn.weit.ktch))
    {
      Change.mean.weight.catch[[s]]=print(fn.plt.mn.ktch.wght(d=Mn.weit.ktch[[s]],NM=capitalize(names(Mn.weit.ktch)[s])))
    }
    figure=ggarrange(plotlist=Change.mean.weight.catch)
    annotate_figure(figure,
                    left = text_grob("Relative mean weight of caught individuals", rot = 90,size=20),
                    bottom = text_grob("Financial year",size=20))
    ggsave(paste(hNdl,'/Outputs/Figure_Changes in mean weight caught individual _TDGDLF.tiff',sep=''),
           width = 10,height = 10,compression = "lzw")
    
  }
  
  clear.log('Change.mean.weight.catch')
  clear.log('Change.mean.length')