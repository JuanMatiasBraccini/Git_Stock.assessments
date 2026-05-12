# 2026 Recover -----------------------------------------------------------------
fn.do.2026.recovery=function(do.fig1=TRUE,do.fig2=TRUE,do.fig3=TRUE,do.fig4=TRUE)
{
  hndl.paper=handl_OneDrive('Scientific manuscripts/Population dynamics/Recovery')
  Paper.species=c("whiskery shark","dusky shark","sandbar shark")
  
  #Figure 1. Management timeline
  if(do.fig1)
  {
    library(patchwork)
    
    source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/Timelines_management.R'))
    p_man.timeline=fun.management.timeline(Management,labl.size=2.5,pt.siz=.9,Right.Margin=0)
    
    #effort 
    built_plot <- ggplot_build(p_man.timeline)
    breaks_numeric <- built_plot$layout$panel_params[[1]]$x$get_breaks()
    breaks_dates <- as.Date(breaks_numeric, origin = "1970-01-01")
    Ef.tdgdlf=Effort.monthly%>%
      mutate(Total=Total/max(Total),
             fishery='TDGDLF',
             year=substr(FINYEAR,1,4))%>%
      dplyr::select(-FINYEAR)
    Ef.nsf=Effort.monthly.north%>%
      rename(Total='Hook days')%>%
      dplyr::select(FINYEAR,Total)%>%
      mutate(Total=Total/max(Total),
             fishery='NSF',
             year=substr(FINYEAR,1,4))%>%
      dplyr::select(-FINYEAR)
    Start.NSF=min(Ef.nsf$year)
    add.ef.yrs=Ef.tdgdlf$year[which(!Ef.tdgdlf$year%in%Ef.nsf$year)]
    ad.dumi=Ef.nsf%>%
      slice(rep(1, length(add.ef.yrs)))%>%
      mutate(Total=0,year=add.ef.yrs)%>%
      filter(year>=Start.NSF)
    Ef.nsf=rbind(Ef.nsf,ad.dumi)%>%arrange(year)
    Ef.dat=rbind(Ef.tdgdlf,Ef.nsf)%>%
      mutate(year=as.POSIXct(paste0("01/01/",year),format="%d/%m/%Y"))
    
    #catch 
    ktch=KtCh%>%
      filter(Name%in%c('dusky shark','sandbar shark','whiskery shark'))%>%
      mutate(fishery=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern shark',
                            ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                     discard.specs),'Southern shark',
                                   ifelse(FishCubeCode%in%c('WRL') & Name%in%WRL.species,'WRL',
                                          'Other'))))
    ktch2=ktch%>%
      group_by(SPECIES,Name,finyear)%>%
      summarise(Tonnes=sum(LIVEWT.c,na.rm=T))%>%
      ungroup()%>%
      group_by(SPECIES,Name)%>%
      mutate(Rel.tonnes=Tonnes/max(Tonnes))%>%
      ungroup()%>%
      mutate(year=as.POSIXct(paste0("01/01/",finyear),format="%d/%m/%Y"),
             Rel.tonnes=SKLR-Rel.tonnes*SKLR,
             Name=capitalize(Name))
    
    #get timeline 
    add.repel.sp.name=FALSE
    LBL.alpha=0.45
    LBL.kl=c('tan4','aquamarine4','navy') 
    Li.wiz=1.25
    p=fun.management.timeline(Management,labl.size=3.15,pt.siz=.9,Right.Margin=120,
                              Start.first.top=TRUE,alpha.decades=0.35,
                              connect.alpha=0.35,connect.kl="black",connect.size=.3,
                              Rev.dec.kl=FALSE,Start.decadal.col='burlywood4')
    #add effort
    add.eFFort=TRUE
    LBL.kl.eff=c('#F8766D','#619CFF')
    names(LBL.kl.eff)=c('NSF','TDGDLF')
    if(add.eFFort)
    {
      LBL.x.pos=ktch2%>%filter(finyear==min(finyear))%>%distinct(year)%>%pull(year)
      LBL.x.pos_NSF=ktch%>%filter(fishery=='Northern shark')%>%filter(finyear==1991)%>%ungroup()%>%
        mutate(year=as.POSIXct(paste0("01/01/",finyear),format="%d/%m/%Y"))%>%distinct(year)%>%pull(year)
      LBL.x.pos_TDGDLF=ktch%>%filter(fishery=='Southern shark' & !FishCubeCode=='Historic')%>%filter(finyear==1975)%>%ungroup()%>%
        mutate(year=as.POSIXct(paste0("01/01/",finyear),format="%d/%m/%Y"))%>%distinct(year)%>%pull(year)
      LBL.y.pos_NSF=min(ktch2$Rel.tonnes)+1 #-26
      LBL.y.pos_TDGDLF=min(ktch2$Rel.tonnes)+1 #-24
      
      Ef.dat2=Ef.dat%>%
        rename(finyear=year)%>%
        group_by(fishery,finyear)%>%
        summarise(Total=sum(Total,na.rm=T))%>%
        ungroup()%>%
        group_by(fishery)%>%
        mutate(Rel.effort=Total/max(Total))%>%
        ungroup()%>%
        mutate(year=as.POSIXct(paste0("01/01/",finyear),format="%d/%m/%Y"),
               Rel.effort=SKLR-Rel.effort*SKLR) 
      NSF.pol=Ef.dat2%>%
        filter(fishery=="NSF")%>%
        dplyr::select(Rel.effort,year)
      NSF.pol=data.frame(year=c(NSF.pol$year,rev(NSF.pol$year)),
                         Rel.effort=c(NSF.pol$Rel.effort,rep(min(NSF.pol$Rel.effort),nrow(NSF.pol))))
      
      TDGDLF.pol=Ef.dat2%>%
        filter(fishery=="TDGDLF")%>%
        dplyr::select(Rel.effort,year)
      TDGDLF.pol=data.frame(year=c(TDGDLF.pol$year,rev(TDGDLF.pol$year)),
                            Rel.effort=c(TDGDLF.pol$Rel.effort,rep(min(NSF.pol$Rel.effort),nrow(TDGDLF.pol))))
      
      
      p=p+
        geom_polygon(data=TDGDLF.pol,aes(year,Rel.effort),
                     alpha=.2,fill=LBL.kl.eff[2])+
        geom_polygon(data=NSF.pol,aes(year,Rel.effort),
                     alpha=.2,fill=LBL.kl.eff[1])+
        #geom_line(data=Ef.dat2%>%filter(fishery=="NSF"),aes(year,Rel.effort),
        #          alpha=LBL.alpha,color=LBL.kl.eff[1],linewidth=1.1,linetype='dashed')+
        #geom_line(data=Ef.dat2%>%filter(fishery=="TDGDLF"),aes(year,Rel.effort),
        #          alpha=LBL.alpha,color=LBL.kl.eff[2],linewidth=1.1,linetype='dashed')+
        geom_text(aes(x=LBL.x.pos_NSF,y=LBL.y.pos_NSF,label="NSF effort"),alpha=LBL.alpha+0.3,color=LBL.kl.eff[1],size=5.5,hjust=0)+
        geom_text(aes(x=LBL.x.pos_TDGDLF,y=LBL.y.pos_TDGDLF,label="TDGDLF effort"),alpha=LBL.alpha+0.3,color=LBL.kl.eff[2],size=5.5,hjust=0)
      
    }
    #add catch
    add.kTch=TRUE
    add.total.ktch.prop=FALSE
    if(add.kTch)
    {
      built_plot <- ggplot_build(p)
      SKLR=built_plot$layout$panel_params[[1]]$y$breaks
      SKLR=subset(SKLR,!is.na(SKLR))
      SKLR=SKLR[1]-7
      p=p+
        geom_line(data=ktch2%>%filter(Name=="Dusky shark"),
                  aes(year,Rel.tonnes),alpha=LBL.alpha,color=LBL.kl[1],linewidth=Li.wiz)+
        geom_line(data=ktch2%>%filter(Name=="Sandbar shark"),
                  aes(year,Rel.tonnes),alpha=LBL.alpha,color=LBL.kl[2],linewidth=Li.wiz)+
        geom_line(data=ktch2%>%filter(Name=="Whiskery shark"),
                  aes(year,Rel.tonnes),alpha=LBL.alpha,color=LBL.kl[3],linewidth=Li.wiz)
      if(add.repel.sp.name)
      {
        p=p+
          geom_text_repel(data=ktch2%>%filter(finyear==1941 & Name=="Dusky shark"),
                          aes(year,Rel.tonnes,label=Name),color=LBL.kl[1],box.padding=10,size=6)+
          geom_text_repel(data=ktch2%>%filter(finyear==1972 & Name=="Sandbar shark"),
                          aes(year,Rel.tonnes,label=Name),color=LBL.kl[2],box.padding=1.8,size=6,hjust=1)+
          geom_text_repel(data=ktch2%>%filter(finyear==1941 & Name=="Whiskery shark"),
                          aes(year,Rel.tonnes,label=Name),color=LBL.kl[3],box.padding=5,size=6,hjust=1)
      }else
      {
        Pos.y=max(ktch2$Rel.tonnes)
        p=p+
          geom_text(data=ktch2%>%filter(finyear==1941 & Name=="Dusky shark"),
                    aes(year,Pos.y-6.5,label=Name),alpha=LBL.alpha+0.1,color=LBL.kl[1],size=6,hjust=0)+
          geom_text(data=ktch2%>%filter(finyear==1941 & Name=="Sandbar shark"),
                    aes(year,Pos.y-10,label=Name),alpha=LBL.alpha+0.1,color=LBL.kl[2],size=6,hjust=0)+
          geom_text(data=ktch2%>%filter(finyear==1941 & Name=="Whiskery shark"),
                    aes(year,Pos.y-3,label=Name),alpha=LBL.alpha+0.1,color=LBL.kl[3],size=6,hjust=0)
      }
      if(add.total.ktch.prop)
      {
        ktch.prop=ktch%>%
          mutate(fishery=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern shark',
                                ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                         discard.specs),'Southern shark',
                                       FishCubeCode)))%>%
          group_by(Name,fishery)%>%
          summarise(Tonnes=sum(LIVEWT.c,na.rm=T))%>%
          arrange(-Tonnes, .by_group = TRUE)%>%
          mutate(Kum=cumsum(Tonnes),
                 Rel=Kum/sum(Tonnes),
                 fishery=ifelse(!Name=='whiskery shark'& Rel<0.9,fishery,
                                ifelse(Name=='whiskery shark'& Rel<0.99,fishery,'other')))%>%
          group_by(Name,fishery)%>%
          summarise(Tonnes=sum(Tonnes,na.rm=T))%>%
          group_by(Name)%>%
          mutate(Prop=Tonnes/sum(Tonnes))%>%
          ungroup()%>%
          mutate(fishery=ifelse(fishery=='Southern shark','TDGDLF',ifelse(fishery=='Northern shark','NSF',fishery)))
        Kl.rec=c("Recreational"="floralwhite");Kl.other=c("other"="#00BA38");Kl.WRL=c("WRL"="black");Kl.taiwan=c('Taiwan'='darkorange')
        p.bar.dusky=ktch.prop%>%
          filter(Name=='dusky shark')%>%
          mutate(fishery=factor(fishery,levels=rev(c('TDGDLF','other','WRL','Recreational'))),
                 LBL=case_when(fishery%in%c('TDGDLF','NSF')~'',
                               !fishery%in%c()~fishery),
                 LBL=capitalize(LBL))%>%
          ggplot(aes(x=1,y=Prop,fill=fishery))+geom_bar(stat='identity',alpha=LBL.alpha,colour = "black")+
          scale_fill_manual(values = c(LBL.kl.eff[2],Kl.rec,Kl.WRL,Kl.other))+
          theme_void()+theme(legend.position = 'top',legend.box.spacing = unit(-15, "pt"),legend.title=element_blank())+
          coord_flip(clip = "off")+
          geom_text(aes(label = LBL),angle = -45,position = position_stack(vjust = 0.5))
        p.bar.sandbar=ktch.prop%>%
          filter(Name=='sandbar shark')%>%
          mutate(fishery=factor(fishery,levels=rev(c('TDGDLF','NSF','other','Taiwan'))),
                 LBL=case_when(fishery%in%c('TDGDLF','NSF','other')~'',
                               !fishery%in%c()~fishery))%>%
          ggplot(aes(x=1,y=Prop,fill=fishery))+geom_bar(stat='identity',alpha=LBL.alpha,colour = "black")+
          scale_fill_manual(values = c(LBL.kl.eff,Kl.taiwan,Kl.other))+
          theme_void()+theme(legend.position = 'top',legend.box.spacing = unit(-15, "pt"),legend.title=element_blank())+
          coord_flip(clip = "off")+
          geom_text(aes(label = LBL),angle = -45,position = position_stack(vjust = 0.5))
        p.bar.whiskery=ktch.prop%>%
          filter(Name=='whiskery shark')%>%
          mutate(fishery=factor(fishery,levels=rev(c('TDGDLF','other'))))%>%
          ggplot(aes(x=1,y=Prop,fill=fishery))+geom_bar(stat='identity',alpha=LBL.alpha,colour = "black")+
          scale_fill_manual(values = c(LBL.kl.eff[2],Kl.other))+
          theme_void()+theme(legend.position = 'top',legend.box.spacing = unit(-15, "pt"),legend.title=element_blank())+
          coord_flip()
        
        y.prop=0.035
        D.y=0.72;S.y=0.63;W.y=0.81
        p=p+
          inset_element(p.bar.dusky, left = 0.04, bottom = D.y, right = 0.4, top = D.y+y.prop)+
          theme(legend.position = 'none')+
          inset_element(p.bar.sandbar, left = 0.04, bottom = S.y, right = 0.4, top = S.y+y.prop)+
          theme(legend.position = 'none')+
          inset_element(p.bar.whiskery, left = 0.04, bottom = W.y, right = 0.4, top = W.y+y.prop)+
          theme(legend.position = 'none')
      }
    }
    #export
    print(p)
    ggsave(paste0(hndl.paper,'/Management & catch.tiff'),width = 11,height = 6, dpi = 300, compression = "lzw")
    
  }
   
  #Figure 2. Catch  
  if(do.fig2)
  {
      SIZ=2
      WID=12
      HEI=10
      t.siz=15
      if(length(Paper.species)>8)
      {
        SIZ=1.5 
        t.siz=13
        WID=14
      }
      if(length(Paper.species)<4) SIZ=3
      Tot.ktch%>%
        filter(Name%in%Paper.species)%>%
        group_by(Type,finyear,Name)%>%
        summarise(Tot=sum(LIVEWT.c,na.rm=T)/unitS)%>%
        filter(!is.na(Tot))%>%
        filter(Tot>0)%>%
        mutate(Name=capitalize(Name))%>%
        ggplot(aes(finyear,Tot,color=Type))+
        geom_point(size=SIZ)+
        geom_line(linetype=2,alpha=0.4)+
        facet_wrap(~Name,scales='free',ncol=1)+
        theme_PA(strx.siz=15,leg.siz=14,axs.t.siz=t.siz,axs.T.siz=18)+
        ylab("Total catch (tonnes)")+xlab("Financial year")+
        theme(legend.position="top",
              legend.title = element_blank(),
              legend.key=element_blank(),
              plot.margin = margin(0.1,.5,0.1,0.1, "cm"))+
        guides(colour = guide_legend(override.aes = list(size=5,linetype = 0)))
    ggsave(paste0(hndl.paper,'/Catch by fleet.tiff'),width = 8,height = 10, dpi = 300, compression = "lzw")
  }
  
  #Figure 3. Biomass and F
  if(do.fig3)
  {
    if(length(Paper.species)>8) InMar=1.25 else InMar=.5
    a=fn.plot.timeseries_combined(this.sp=Paper.species,
                                  d=Age.based$SS,
                                  YLAB="Relative biomass",
                                  Type="Depletion",
                                  InnerMargin=InMar,
                                  RefPoint=Ref.points,
                                  Kach=Age.based$SS$rel.biom)
    b=fn.plot.timeseries_combined(this.sp=Paper.species,
                                  d=Age.based$SS,
                                  YLAB=expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),
                                  Type="F.series",
                                  InnerMargin=InMar,
                                  RefPoint=NULL,
                                  Kach=Age.based$SS$rel.biom)
    ggarrange(a,b,widths = c(0.55, 0.45))
    ggsave(paste0(hndl.paper,'/Relative biomass and F.tiff'),width = 10,height = 9, dpi = 300, compression = "lzw")
  }
  
  #Figure 4. Risk
  if(do.fig4)
  {
    Final.risk.table=fn.risk.figure(d=Weighted.overall.risk_future%>%
                     filter(tolower(Species)%in%Paper.species),
                   Risk.colors=RiskColors,
                   out.plot=TRUE)
    ggsave(paste0(hndl.paper,'/Risk_final.tiff'),width = 10,height = 11, dpi = 300, compression = "lzw")
    write.csv(Final.risk.table,paste0(hndl.paper,'/Risk_final.csv'),row.names=F)
  }
}  

# Sawfish -----------------------------------------------------------------
fn.do.sawfish=function()
{
  hNdl.sawfish=handl_OneDrive(paste('Analyses/Population dynamics/Sawfishes/',Year.of.assessment,sep=''))
  if(!dir.exists(hNdl.sawfish))dir.create(hNdl.sawfish)
  
  #1. Plot catches
  Sawfish.ktch=Get.ktch$Total.method%>%
    filter(SPECIES%in%25000:25020)%>%
    filter(finyear>=1960)%>%
    mutate(Name=capitalize(Name))
  
  xx=Sawfish.ktch%>%
    data.frame%>%
    distinct(Name,SPECIES)
  assessed.sawfish=xx%>%pull(Name)
  names(assessed.sawfish)=xx%>%pull(SPECIES)
  n.sawfish=length(assessed.sawfish)
  
  Sawfish.ktch%>%
    filter(LIVEWT.c>0)%>%
    mutate(Fishery.gear=paste(Gear,FishCubeCode,sep='-'))%>%
    ggplot(aes(finyear,LIVEWT.c,colour=Fishery.gear))+
    geom_point()+geom_line(alpha=0.25)+
    facet_wrap(~Name,nrow=3,scales='free')+
    theme_PA(strx.siz=14,leg.siz=11,axs.t.siz=13,axs.T.siz=16)+
    ylab("Catch (tonnes)")+xlab("Financial year")+
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=5,linetype = 0)))
  ggsave(paste(hNdl.sawfish,'Annual_ktch_by_species.tiff',sep='/'), width = 10,height = 10, dpi = 300, compression = "lzw")
  
  for(i in 1:n.sawfish)  #missing
  {
    #DBSRA
    yrs=Catch_only_sawfish$DBSRA[[i]][['S1']]$output$Years
    Bmsy=apply(Catch_only_sawfish$DBSRA[[i]][['S1']]$output$B.Bmsy,2,median,na.rm=T)
    Fmsy=apply(Catch_only_sawfish$DBSRA[[i]][['S1']]$output$F.Fmsy,2,median,na.rm=T)
    p.DBSRA=kobePlot(f.traj=Fmsy[1:length(yrs)],
                     b.traj=Bmsy[1:length(yrs)],
                     Years=yrs,
                     Titl=paste("DBSRA",names(Catch_only_sawfish$DBSRA)[i],sep='-'))
    rm(yrs,Fmsy,Bmsy)
    
    #CMSY
    if(CMSY.method=="Froese")
    {
      yrs=Catch_only_sawfish$CMSY[[i]][['S1']]$output$ref_ts$year
      Bmsy=Catch_only_sawfish$CMSY[[i]][['S1']]$output$ref_ts$bbmsy
      Fmsy=Catch_only_sawfish$CMSY[[i]][['S1']]$output$ref_ts$ffmsy
      
    }
    if(CMSY.method=="Haddon")
    {
      yrs=Catch_only_sawfish$CMSY[[i]][['S1']]$output$Years
      Bmsy=apply(Catch_only_sawfish$CMSY[[i]][['S1']]$output$B.Bmsy,2,median,na.rm=T)[1:length(yrs)]
      Fmsy=apply(Catch_only_sawfish$CMSY[[i]][['S1']]$output$F.Fmsy,2,median,na.rm=T)[1:length(yrs)]
    }
    
    p.CMSY=kobePlot(f.traj=Fmsy,
                    b.traj=Bmsy,
                    Years=yrs,
                    Titl=paste("CMSY",names(Catch_only_sawfish$DBSRA)[i],sep='-'))
    rm(yrs,Fmsy,Bmsy)
    
    #JABBA
    p.JABBA=with(Catch_only_sawfish$JABBA[[i]][['S1']]$output,
                 {
                   kobePlot(f.traj=timeseries[, , "FFmsy"][,"mu"],
                            b.traj=timeseries[, , "BBmsy"][,"mu"],
                            Years=yr,
                            Titl=paste("JABBA",names(JABBA.sawfish)[i],sep='-'),
                            Probs=data.frame(x=kobe$stock,
                                             y=kobe$harvest))
                 })
    
    figure <- ggarrange(plotlist=list(p.DBSRA+rremove("axis.title"),
                                      p.CMSY+rremove("axis.title"),
                                      p.JABBA+rremove("axis.title")),
                        ncol=1,nrow=3,common.legend = FALSE)
    
    annotate_figure(figure,
                    bottom = text_grob(expression(B/~B[MSY]), size=16),
                    left = text_grob(expression(F/~F[MSY]), rot = 90,size=16))
    
    ggsave(paste(hNdl.sawfish,paste('Kobe_',names(Catch_only_sawfish$DBSRA)[i],'.tiff',sep=''),sep='/'),
           width = 8,height = 14, dpi = 300, compression = "lzw")
    
    
  }
  
}