#---Spatio-temporal catch and effort. Reported TDGLF and NSF ----   
#note: bubble size is proportion of blocks fished out of maximum number of blocks fished for each species

#get reported catch
dis.sp=All.species.names%>%filter(SNAME%in%unlist(Lista.sp.outputs))%>%pull(SPECIES)
dis.sp=c(dis.sp,19000)
Southern=fn.in(NM='Data.monthly.csv')%>%
  filter(!Shark.fishery=='non.shark.fishery' & SPECIES%in%dis.sp)%>%
  filter(!is.na(BLOCKX))%>%
  group_by(FINYEAR,BLOCKX,SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
Northern=fn.in(NM='Data.monthly.NSF.csv')%>%
  filter(!Shark.fishery=='non.shark.fishery' & SPECIES%in%dis.sp)%>%
  filter(!is.na(BLOCKX))%>%
  group_by(FINYEAR,BLOCKX,SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
Spatio.temp.dat=rbind(Southern%>%
                        mutate(Fishery='Southern'),
                      Northern%>%
                        mutate(Fishery='Northern'))
rm(Southern,Northern)

Relevant.yrs=sort(unique(KtCh.method$FINYEAR))

#spatial effort
FINYrS=sort(unique(c(as.character(unique(Effort.monthly_blocks$FINYEAR))),sort(as.character(unique(Effort.daily_blocks$finyear)))))
grouping=5
FINYrS.gp=seq(1,length(FINYrS),by=grouping)
FINYrS.gped=vector('list',length(FINYrS.gp))
for(f in 1:length(FINYrS.gped))
{
  if(f==length(FINYrS.gped))
  {
    FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:length(FINYrS)]
    if(length(FINYrS.gped[[f]])==1)
    {
      names(FINYrS.gped)[f]=FINYrS.gped[[f]][1]
    }else
    {
      names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
    }
    
  }else
  {
    FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:(FINYrS.gp[f+1]-1)]
    names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
  }
}
FINYrS.gped=FINYrS.gped[which(unlist(lapply(FINYrS.gped,function(x) length(x)==5)))]
FINYrS.gped=do.call(cbind,FINYrS.gped)%>%
  data.frame%>%
  gather('Rango','FINYEAR')%>%
  mutate(Rango=substr(Rango,2,50),
         Rango=sub(".to.", " to ", Rango),
         Rango=str_replace_all(Rango, c("\\."), "-"))
daily.years=unique(Effort.daily_blocks$finyear)
daily.years=subset(daily.years,!daily.years=="2005-06")
Spatial.effort_monthly=Effort.monthly_blocks%>%
  filter(!Shark.fishery=='non.shark.fishery' &
           FINYEAR%in%Relevant.yrs&
           !FINYEAR%in%daily.years)%>%
  mutate(Effort=Km.Gillnet.Days.c,
         LAT1=-floor(abs(LAT)),
         LONG1=floor(LONG))%>%
  filter(!is.na(Effort) | !is.na(LAT))%>%
  filter(METHOD=='GN')%>%
  group_by(Same.return,LAT1,LONG1,FINYEAR)%>%
  summarise(Effort=max(Effort))%>%
  ungroup()%>%
  group_by(LAT1,LONG1,FINYEAR)%>%
  summarise(Effort=sum(Effort))
Spatial.effort_daily=Effort.daily_blocks%>%
  filter(!Shark.fishery=='non.shark.fishery' &
           finyear%in%Relevant.yrs )%>%
  mutate(Effort=Km.Gillnet.Days.c,
         LAT1=-floor(abs(LAT)),
         LONG1=floor(LONG))%>%
  filter(!is.na(Effort) | !is.na(LAT))%>%
  filter(method=='GN')%>%
  group_by(Same.return.SNo,LAT1,LONG1,finyear)%>%
  summarise(Effort=max(Effort))%>%
  ungroup()%>%
  rename(FINYEAR=finyear)%>%
  group_by(LAT1,LONG1,FINYEAR)%>%
  summarise(Effort=sum(Effort))
Spatial.effort_monthly.north=Effort.monthly.north_blocks%>%
  filter(FINYEAR%in%Relevant.yrs &
           !FINYEAR%in%unique(Effort.daily.north_blocks$finyear))%>%
  mutate(Effort=hook.days,
         LAT1=-floor(abs(as.numeric(substr(BLOCKX,1,2)))),
         LONG1=floor(100+as.numeric(substr(BLOCKX,3,4))),
         LAT1=ifelse(BLOCKX%in%c(96021),-25,LAT1),
         LONG1=ifelse(BLOCKX%in%c(96021),113,LONG1))%>%
  filter(!is.na(Effort) | !is.na(LAT1))%>%
  filter(METHOD=='LL')%>%
  group_by(Same.return,LAT1,LONG1,FINYEAR)%>%
  summarise(Effort=max(Effort))%>%
  ungroup()%>%
  group_by(LAT1,LONG1,FINYEAR)%>%
  summarise(Effort=sum(Effort))
Spatial.effort_daily.north=Effort.daily.north_blocks%>%
  rename(METHOD=method,
         FINYEAR=finyear,
         BLOCKX=blockx)%>%
  filter(FINYEAR%in%Relevant.yrs)%>%
  mutate(Effort=hook.days,
         LAT1=-floor(abs(as.numeric(substr(BLOCKX,1,2)))),
         LONG1=floor(100+as.numeric(substr(BLOCKX,3,4))),
         LAT1=ifelse(BLOCKX%in%c(96021),-25,LAT1),
         LONG1=ifelse(BLOCKX%in%c(96021),113,LONG1))%>%
  filter(!is.na(Effort) | !is.na(LAT1))%>%
  filter(METHOD=='LL')%>%
  group_by(Same.return.SNo,LAT1,LONG1,FINYEAR)%>%
  summarise(Effort=max(Effort))%>%
  ungroup()%>%
  group_by(LAT1,LONG1,FINYEAR)%>%
  summarise(Effort=sum(Effort))
Spatial.effort=rbind(Spatial.effort_monthly%>%mutate(Fishery='Southern'),
                     Spatial.effort_daily%>%mutate(Fishery='Southern'),
                     Spatial.effort_monthly.north%>%mutate(Fishery='Northern'),
                     Spatial.effort_daily.north%>%mutate(Fishery='Northern'))%>%
  filter(!is.na(Effort))%>%
  left_join(FINYrS.gped,by='FINYEAR')%>%
  group_by(LAT1,LONG1,Rango,Fishery)%>%
  summarise(Effort=sum(Effort))%>%
  group_by(Fishery)%>%
  mutate(Rel.effort=Effort/max(Effort))%>%
  ungroup()

#number of fished blocks
Effort.monthly_blocks=Effort.monthly_blocks%>%
  filter(FINYEAR%in%Relevant.yrs)%>%
  count(FINYEAR,BLOCKX)%>%
  group_by(FINYEAR,BLOCKX)%>%
  mutate(n=ifelse(n>0,1,0))
Effort.daily_blocks=Effort.daily_blocks%>%
  rename(FINYEAR=finyear,
         BLOCKX=blockx)%>%
  filter(FINYEAR%in%Relevant.yrs)%>%
  count(FINYEAR,BLOCKX)%>%
  group_by(FINYEAR,BLOCKX)%>%
  mutate(n=ifelse(n>0,1,0))
Effort_blocks=rbind(Effort.monthly_blocks,Effort.daily_blocks)%>%
  count(FINYEAR,BLOCKX)%>%
  group_by(FINYEAR,BLOCKX)%>%
  mutate(n=ifelse(n>0,1,0))%>%
  group_by(FINYEAR)%>%
  summarise(Tot=sum(n))%>%
  data.frame

Effort.monthly.north_blocks=Effort.monthly.north_blocks%>%
  filter(FINYEAR%in%Relevant.yrs)%>%
  count(FINYEAR,BLOCKX)%>%
  group_by(FINYEAR,BLOCKX)%>%
  mutate(n=ifelse(n>0,1,0))
Effort.daily.north_blocks=Effort.daily.north_blocks%>%
  rename(FINYEAR=finyear,
         BLOCKX=blockx)%>%
  filter(FINYEAR%in%Relevant.yrs)%>%
  count(FINYEAR,BLOCKX)%>%
  group_by(FINYEAR,BLOCKX)%>%
  mutate(n=ifelse(n>0,1,0))
Effort.north_blocks=rbind(Effort.monthly.north_blocks,Effort.daily.north_blocks)%>%
  count(FINYEAR,BLOCKX)%>%
  group_by(FINYEAR,BLOCKX)%>%
  mutate(n=ifelse(n>0,1,0))%>%
  group_by(FINYEAR)%>%
  summarise(Tot=sum(n))%>%
  data.frame

NSF.no.fishing=data.frame(FINYEAR=Effort_blocks$FINYEAR[which(!Effort_blocks$FINYEAR%in%Effort.north_blocks$FINYEAR)],
                          Tot=0)
Effort_blocks.all=rbind(Effort_blocks%>%mutate(Fishery='Southern'),
                        rbind(Effort.north_blocks,NSF.no.fishing)%>%
                          arrange(FINYEAR)%>%
                          mutate(Fishery='Northern'))


#plot
fn.spatio.temp.catch.dist=function(d,Snames,prop.by='fishery',show.prop='blks.with.ktch.over.fished')
{
  #some manipulations
  this.sp=All.species.names%>%filter(SNAME%in%Snames)
  
  d1=d%>%
    filter(SPECIES%in%this.sp$SPECIES)%>%
    left_join(this.sp,by='SPECIES')%>%
    ungroup()
  Unik.yr=d1%>%
    distinct(FINYEAR)%>%
    arrange(FINYEAR)%>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))
  Unik.yr$id.x <- 1:nrow(Unik.yr) 
  
  Unik.sp=d1%>%
    distinct(SNAME)%>%
    arrange(SNAME)
  Unik.sp$id.y <- 1:nrow(Unik.sp) 
  
  Xlab=Unik.yr$year
  names(Xlab)=Unik.yr$id.x
  Xlab=Xlab[seq(1,length(Xlab),10)]
  
  Ylab=capitalize(Unik.sp$SNAME)
  names(Ylab)=Unik.sp$id.y
  
  Dummy=Effort_blocks.all%>%
    left_join(Unik.yr,by='FINYEAR')
  
  #spatio-temporal effort distribution
  p1=Spatial.effort%>%
    ggplot(aes(LONG1,LAT1))+
    geom_raster(aes(fill = Rel.effort))+
    facet_wrap(~Rango,ncol=2)+
    scale_fill_gradientn(colours=c('ivory2','gold',"red2","darkred"))+
    theme_PA(leg.siz=12.5,axs.t.siz=13,axs.T.siz=19,Sbt.siz=14,strx.siz=14)+
    theme(
      legend.title=element_blank(),
      legend.key.size = unit(1, 'cm'),
      legend.position = c(0.8, 0.075))+
    ylab('Latitude (S)')+xlab('Longitude (E)')
  
  
  #annual blocks with catch over blocks fished
  Max.ktch=d1%>%
    group_by(SNAME)%>%
    summarise(Ktch.mx=sum(LIVEWT.c,na.rm=T))%>%
    ungroup()
  
  ktch.yr.sp=d1%>%
    group_by(SNAME,FINYEAR)%>%
    summarise(Ktch=sum(LIVEWT.c,na.rm=T))%>%
    ungroup()%>%
    left_join(Max.ktch,by='SNAME')%>%
    mutate(prop.ktch=Ktch/Ktch.mx)%>%
    dplyr::select(-c(Ktch,Ktch.mx))
  
  if(show.prop=='fished.blks.with.ktch')
  {
    d2=d1%>%
      count(FINYEAR,SNAME,BLOCKX,Fishery)%>%
      group_by(FINYEAR,SNAME,Fishery)%>%
      mutate(n=ifelse(n>0,1,0))%>%
      group_by(FINYEAR,SNAME,Fishery)%>%
      summarise(n=sum(n,na.rm=T))%>%
      ungroup
    if(prop.by=='fishery')
    {
      Max.blk.yr=d1%>%
        distinct(BLOCKX,Fishery)%>%
        count(BLOCKX,Fishery)%>%
        group_by(Fishery)%>%
        mutate(n=ifelse(n>0,1,0))%>%
        group_by(Fishery)%>%
        summarise(Max.n=sum(n,na.rm=T))%>%
        ungroup
      
      d2=d2%>%
        left_join(Max.blk.yr,by=c('Fishery'))%>%
        mutate(prop=n/Max.n)%>%
        ungroup%>%
        data.frame
    }
    if(prop.by=='year.fishery')
    {
      Max.blk.yr=d1%>%
        count(FINYEAR,BLOCKX,Fishery)%>%
        group_by(FINYEAR,Fishery)%>%
        mutate(n=ifelse(n>0,1,0))%>%
        group_by(FINYEAR,Fishery)%>%
        summarise(Max.n=sum(n,na.rm=T))%>%
        ungroup
      
      d2=d2%>%
        left_join(Max.blk.yr,by=c('FINYEAR','Fishery'))%>%
        mutate(prop=n/Max.n)%>%
        ungroup%>%
        data.frame
    }
    
    d2=d2%>%
      left_join(ktch.yr.sp,by=c('SNAME','FINYEAR'))%>%
      filter(!(SNAME%in%c('gummy shark','whiskery shark') & Fishery=='Northern' ))
    
  }
  if(show.prop=='blks.with.ktch.over.fished')
  {
    N.fished.blocks.yr=d1%>%
      distinct(FINYEAR,BLOCKX,Fishery)%>%
      count(FINYEAR,Fishery)%>%
      rename(Fished.blks=n)
    
    N.blocks.with.ktch.yr=d1%>%
      count(FINYEAR,SNAME,BLOCKX,Fishery)%>%
      group_by(FINYEAR,SNAME,Fishery)%>%
      mutate(n=ifelse(n>0,1,0))%>%
      group_by(FINYEAR,SNAME,Fishery)%>%
      summarise(n=sum(n,na.rm=T))%>%
      ungroup
    
    d2=left_join(N.blocks.with.ktch.yr,
                 N.fished.blocks.yr,by=c('FINYEAR','Fishery'))%>%
      mutate(Prop.ktch.over.fishd=n/Fished.blks)%>%
      group_by(SNAME,Fishery)%>%
      mutate(Max=max(Prop.ktch.over.fishd))%>%
      ungroup()%>%
      mutate(prop=Prop.ktch.over.fishd/Max)
  }
  d2=d2%>%
    filter(Fishery=='Southern')  #not applicable to 'Northern' as ceased operations in 2008-09
  
  coeff=max(Dummy$Tot)
  p2=d2%>%
    mutate(yr=as.numeric(substr(FINYEAR,1,4)),
           SNAME=capitalize(SNAME))%>%
    ggplot(aes(yr,prop,colour=Fishery))+
    geom_point(size=2,show.legend = FALSE)+
    #geom_line(linetype=1,size=1)+
    facet_wrap(~SNAME,ncol=2)+
    ylab('Proportion of fished blocks with catch')+
    xlab('Financial year')+
    theme_PA(axs.t.siz=13,axs.T.siz=19,Sbt.siz=14,leg.siz=20,strx.siz=16)+
    ylim(0,1)+
    theme(legend.position = 'none',
          legend.title=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=2)))+
    geom_line(data=Dummy%>%filter(Fishery=='Southern'), 
              aes(x=year,y=Tot/coeff),size=1.25,color='black',alpha=0.25)+
    # geom_line(data=Dummy%>%filter(Fishery=='Northern'),
    #            aes(x=year,y=Tot/coeff),size=1,color='black',alpha=0.5,linetype = 'longdash')+
    scale_y_continuous(sec.axis = sec_axis(~.*coeff, name="Total number of fished blocks"))+
    theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
    theme(axis.title.y.left = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
  
  #combine figures
  figure=ggarrange(plotlist=list(p1,p2),ncol=2,heights=c(1,1))+
    theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
  print(figure)
  
  #annual blocks with catch over blocks ever catching the species 
  Max.blk.sp=d1%>%
    count(SNAME,BLOCKX)%>%
    group_by(SNAME)%>%
    mutate(n=ifelse(n>0,1,0))%>%
    group_by(SNAME)%>%
    summarise(Max.n=sum(n,na.rm=T))%>%
    ungroup
  
  d1=d1%>%
    count(FINYEAR,SNAME,BLOCKX)%>%
    group_by(FINYEAR,SNAME)%>%
    mutate(n=ifelse(n>0,1,0))%>%
    group_by(FINYEAR,SNAME)%>%
    summarise(n=sum(n,na.rm=T))%>%
    ungroup%>%
    left_join(Max.blk.sp,by='SNAME')%>%
    mutate(prop=n/Max.n)%>%
    ungroup%>%
    data.frame%>%
    left_join(Unik.yr,by='FINYEAR')%>%
    left_join(Unik.sp,by='SNAME')
  
  add.to.ylab=d1%>%distinct(SNAME,Max.n)%>%mutate(SNAME=capitalize(SNAME))
  add.to.ylab=add.to.ylab%>%
    left_join(data.frame(SNAME=Ylab,index=as.numeric(names(Ylab))),
              by='SNAME')%>%
    arrange(index)%>%
    mutate(LBL=paste(SNAME,' (n=',Max.n,')',sep=''))
  
  Ylab=add.to.ylab$LBL
  names(Ylab)=add.to.ylab$index
  
  do.this=FALSE
  if(do.this)
  {
    coeff=max(Dummy$Tot)/max(Unik.sp$id.y)
    
    p3=d1%>%
      ggplot(aes(x=year,y=id.y)) +
      geom_point(aes(size=prop),colour='orange')+
      ylab('')+xlab('')+
      theme_PA(axs.t.siz=13,axs.T.siz=16,Sbt.siz=14)+
      theme(legend.position='none',
            plot.subtitle = element_text(vjust = 1))+
      geom_line(data=Dummy%>%filter(Fishery=='Southern'), aes(x=year,y=Tot/coeff),size=1.5,
                color='blue',alpha=0.6,linetype = 1)+
      geom_line(data=Dummy%>%filter(Fishery=='Northern'), aes(x=year,y=Tot/coeff),size=1.5,
                color='red',alpha=0.6,linetype = 1)+
      scale_x_continuous(labels=Xlab,breaks=as.numeric(names(Xlab)))+
      scale_y_continuous(labels=Ylab,breaks=as.numeric(names(Ylab)),
                         sec.axis = sec_axis(~.*coeff, name="Number of blocks fished"))+
      labs(subtitle='Blocks with landings over the maximum number of blocks with landings for each species')
    # scale_colour_gradient(low = "darkred", high = "darkgoldenrod1", na.value = NA)
    
  }
  
  # DD1=d1%>%
  #   dplyr::select(FINYEAR,prop,SNAME)%>%
  #   spread(FINYEAR,prop,fill=0)
  # DD=as.matrix(DD1%>%dplyr::select(-SNAME))
  # rownames(DD)=DD1$SNAME
  return(list(Fished.blks.Fishery=Dummy,
              prop.ktch_over.fished.bocks=d2))
}
Store.spatial.temporal.ktch=Lista.sp.outputs[-match('additional.sp',names(Lista.sp.outputs))]
for(l in 1:length(Store.spatial.temporal.ktch))
{
  wiz = 15
  heiz = 12
  get.dis.sp=Lista.sp.outputs[[l]]
  if(names(Lista.sp.outputs)[l]=="Other.sp")
  {
    get.dis.sp=c(get.dis.sp,'hammerheads')
    wiz = 12
    heiz = 14
  }
  Store.spatial.temporal.ktch[[l]]=fn.spatio.temp.catch.dist(d=Spatio.temp.dat,
                                                             Snames=get.dis.sp)
  ggsave(paste(Rar.path,'/Spatio.temporal.catch_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = wiz,height = heiz,compression = "lzw")
}
clear.log('fn.spatio.temp.catch.dist')
clear.log('Spatio.temp.dat')
clear.log('Effort.daily_blocks')
clear.log('Effort.monthly_blocks')

