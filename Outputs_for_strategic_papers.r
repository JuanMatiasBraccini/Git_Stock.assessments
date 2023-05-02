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