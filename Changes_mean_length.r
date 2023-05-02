  #23.1   TDGDLF
  Change.mean.length=vector('list',N.sp)
  names(Change.mean.length)=names(Species.data)
  fun.change.mean.len=function(d,NM,toMatch,min.annual.obs,XLIM,len.mat,max.len)
  {
    d.list=d[grep(paste(toMatch,collapse="|"),names(d))]
    if(sum(grepl('Table',names(d.list)))>0)d.list=d.list[-grep('Table',names(d.list))]
    
    if(length(d.list)>0)
    {
      my_formula = y ~ x
      
      #by mesh
      d.list=do.call(rbind,d.list)
      d.list=d.list%>%
        mutate(mesh=case_when(grepl("6.5.inch",rownames(d.list))~6.5,
                              grepl("7.inch",rownames(d.list))~7))
      N.min=d.list%>%
        group_by(FINYEAR,mesh)%>%  
        tally()%>%
        filter(n>=min.annual.obs)%>%
        mutate(Keep="YES")
      
      d.list=d.list%>%
        left_join(N.min,by=c('FINYEAR','mesh'))%>%
        filter(Keep=="YES")
      if(nrow(d.list)>0 & length(unique(d.list$FINYEAR))>2)
      {
        p=d.list%>%
          mutate(Finyear=as.numeric(substr(FINYEAR,1,4)),
                 Finyear.d=factor(Finyear,levels=sort(unique(Finyear))),
                 Mesh=paste(mesh,'inch'))%>%      
          ggplot(aes(x = Finyear, y = FL)) +
          geom_violin(aes(fill = Finyear, group = Finyear.d), alpha = 0.2) +
          facet_wrap(~Mesh,scales='free')+
          stat_summary(fun = "mean",geom = "point",color = "red",size=2)+
          geom_smooth(color = rgb(1,.1,.1,alpha=0.2), formula = my_formula, method = 'lm',se=TRUE)+ 
          stat_poly_eq(aes(label = paste("atop(", stat(eq.label),  ",", 
                                         paste(stat(adj.rr.label),stat(p.value.label), sep = "*\", \"*"), ")")),
                       formula = my_formula, parse = TRUE,
                       label.y = "top", label.x = "right", size = 4) +
          xlab("")+ylab("")+
          theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=16)+
          theme(legend.position = "none",
                plot.title =element_text(size=17))+
          labs(title=NM)+
          xlim(XLIM)+
          ylim(0,max.len)+
          geom_hline(yintercept = max.len,colour='orange',size=1.05)+
          geom_hline(yintercept = len.mat,colour='orange',linetype = "dashed",size=1.05)
        #ylim(quantile(d.list$FL,na.rm=T,probs = 0.01),quantile(d.list$FL,na.rm=T,probs = 0.975))  
        return(p)
      }
    }
  }
  for(s in 1:N.sp) 
  {
    print(paste("Change in mean length ","--",names(Species.data)[s]))
    dummy=print(fun.change.mean.len(d=Species.data[[s]],
                                    NM=capitalize(names(Species.data)[s]),
                                    toMatch=c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2"),
                                    min.annual.obs=Min.annual.obs.ktch,
                                    XLIM=c(1990,as.numeric(substr(Last.yr.ktch,1,4))),
                                    len.mat=(List.sp[[s]]$TL.50.mat-List.sp[[s]]$b_FL.to.TL)/List.sp[[s]]$a_FL.to.TL,  #total length to fl
                                    max.len=(List.sp[[s]]$TLmax-List.sp[[s]]$b_FL.to.TL)/List.sp[[s]]$a_FL.to.TL))
    if(!is.null(dummy)) Change.mean.length[[s]]=dummy
    rm(dummy)
  }
  
  
  #23.2   Survey
  do.survey.size=FALSE
  if(do.survey.size)
  {
    Survey.mean.size=vector('list',N.sp)
    names(Survey.mean.size)=Keep.species
    for(l in 1:N.sp)
    {
      dummy=list()
      if('Srvy.FixSt_size'%in%names(Species.data[[l]]))
      {
        d=Species.data[[l]]$Srvy.FixSt_size%>%
          mutate(Finyear=paste(yr-1,substr(yr,3,4),sep='-'),
                 yr.f=as.numeric(substr(Finyear,1,4)))%>%
          rename(Mean=MeAn,
                 LOW.CI=LowCI,
                 UP.CI=UppCI)%>%
          filter(yr.f<=Last.yr.ktch.numeric)
        if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$Survey=d
        rm(d)
      }
      if(length(dummy)>0) Survey.mean.size[[l]]=dummy
    }
    fun.change.mean.len.survey=function(d,NM,Ymin='data',Xlim)
    {
      if(!is.null(d))
      {
        if(Ymin=='data') YLIM.min=min(d$Survey$LOW.CI) else YLIM.min=0
        my_formula = y ~ x
        p=d$Survey%>%
          mutate(Mesh='Survey')%>%
          ggplot(aes(yr.f,Mean))+
          geom_point(aes(color=yr.f),size=2,show.legend=F)+
          geom_errorbar(aes(ymin=UP.CI, ymax=LOW.CI,color=yr.f), width=.5)+
          facet_wrap(~Mesh,scales='free')+
          geom_smooth(color = "black",alpha=0.3, formula = my_formula, method = 'lm',se=FALSE)+ 
          stat_poly_eq(aes(label = paste("atop(", stat(eq.label),  ",", 
                                         paste(stat(adj.rr.label),stat(p.value.label), sep = "*\", \"*"), ")")),
                       formula = my_formula, parse = TRUE,
                       label.y = "bottom", label.x = "right", size = 4) +
          xlab("")+ylab("")+
          theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=16)+
          theme(legend.position = "none",
                plot.title =element_text(size=17))+
          labs(title=NM)+
          xlim(Xlim)+
          ylim(YLIM.min,max(d$Survey$UP.CI))
        return(p)
      }
    } 
    Change.mean.length_Survey=vector('list',N.sp)
    names(Change.mean.length_Survey)=Keep.species
    for(l in 1:N.sp) 
    {
      print(paste("Change in mean length ","--",names(Survey.mean.size)[l]))
      dummy=print(fun.change.mean.len.survey(d=Survey.mean.size[[l]],
                                             NM=capitalize(names(Survey.mean.size)[l]),
                                             Ymin='data',
                                             Xlim=c(2000,Last.yr.ktch.numeric+1)))
      if(!is.null(dummy)) Change.mean.length_Survey[[l]]=dummy
    }
    for(l in 1:length(Store.spatial.temporal.ktch))
    {
      figure=ggarrange(plotlist=fun.find.in.list(x=Change.mean.length_Survey[Lista.sp.outputs[[l]]]),
                       ncol=1)+
        theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
      annotate_figure(figure,
                      left = text_grob("Relative size (95% CI)", rot = 90,size=20,vjust=1),
                      bottom = text_grob("Financial year",size=20,vjust=-1))
      
      ggsave(paste(Rar.path,'/Changes in observed mean length_Survey_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
             width = 10,height = 10,compression = "lzw")
    }
  }
