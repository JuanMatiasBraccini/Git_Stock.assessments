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
Northern=fn.in(NM='Data.monthly.north.csv')%>%
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

#Keep only relevant years
ol.iers=sort(unique(Spatio.temp.dat$FINYEAR))
use.dis.yrs=ol.iers[1:match(Last.yr.ktch,ol.iers)]
Spatio.temp.dat=Spatio.temp.dat%>%filter(FINYEAR%in%use.dis.yrs)


#spatial effort
Effort.monthly_blocks=fn.in("Effort.monthly.csv")
Effort.daily_blocks=fn.in("Effort.daily.csv")
Effort.monthly.north_blocks=fn.in("Effort.monthly.north.csv")
Effort.daily.north_blocks=fn.in("Effort.daily.north.csv") 

Effort.monthly_blocks=Effort.monthly_blocks%>%filter(FINYEAR%in%use.dis.yrs) 
Effort.monthly.north_blocks=Effort.monthly.north_blocks%>%filter(FINYEAR%in%use.dis.yrs) 
Effort.daily_blocks=Effort.daily_blocks%>%filter(finyear%in%use.dis.yrs) 
Effort.daily.north_blocks=Effort.daily.north_blocks%>%filter(finyear%in%use.dis.yrs) 

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
remove.incomplete.years=FALSE
if(remove.incomplete.years)FINYrS.gped=FINYrS.gped[which(unlist(lapply(FINYrS.gped,function(x) length(x)==5)))]
suppressWarnings({FINYrS.gped=do.call(cbind,FINYrS.gped)%>%
  data.frame%>%
  gather('Rango','FINYEAR')%>%
  mutate(Rango=substr(Rango,2,50),
         Rango=sub(".to.", " to ", Rango),
         Rango=str_replace_all(Rango, c("\\."), "-"))%>%
  distinct(Rango,FINYEAR)})
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

BLOCKX_lat_long=rbind(Effort.monthly_blocks%>%
                        distinct(LAT,LONG,BLOCKX)%>%
                        mutate(LAT=-(abs(LAT)),
                               LONG=(LONG)),
                      Effort.monthly.north_blocks%>%
                        mutate(LAT=-floor(abs(as.numeric(substr(BLOCKX,1,2)))),
                               LONG=floor(100+as.numeric(substr(BLOCKX,3,4))),
                               LAT=ifelse(BLOCKX%in%c(96021),-25,LAT),
                               LONG=ifelse(BLOCKX%in%c(96021),113,LONG))%>%
                        distinct(LAT,LONG,BLOCKX)%>%
                        mutate(LAT=-(abs(LAT)),
                               LONG=(LONG)),
                      Effort.daily_blocks%>%
                        rename(BLOCKX=blockx)%>%
                        distinct(LAT,LONG,BLOCKX)%>%
                        mutate(LAT=-floor(abs(LAT)),
                               LONG=floor(LONG))%>%
                        distinct(BLOCKX,.keep_all = T),
                      Effort.daily.north_blocks%>%
                        rename(BLOCKX=blockx)%>%
                        mutate(LAT=-floor(abs(as.numeric(substr(BLOCKX,1,2)))),
                               LONG=floor(100+as.numeric(substr(BLOCKX,3,4))),
                               LAT=ifelse(BLOCKX%in%c(96021),-25,LAT),
                               LONG=ifelse(BLOCKX%in%c(96021),113,LONG))%>%
                        distinct(LAT,LONG,BLOCKX)%>%
                        mutate(LAT=-(abs(LAT)),
                               LONG=(LONG))%>%
                        distinct(BLOCKX,.keep_all = T))%>%
  filter(!is.na(LAT))%>%
  distinct(BLOCKX,.keep_all = T)

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
                          mutate(Fishery='Northern'))%>%
                  filter(FINYEAR%in%use.dis.yrs)


#get species CSIRO distribution
library(rgdal)
library(sf)
library(ozmaps)
source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/CSIRO distribution maps.r'))  
hendl=handl_OneDrive('Analyses/Mapping/Map reported species distribution/')
CAAB=read.csv(paste(hendl,'caab_dump_latest.csv',sep=''))%>%
        filter(CLASS=='Elasmobranchii')%>%
        mutate(COMMON_NAME=ifelse(COMMON_NAME=='Fossil Shark','Snaggletooth',COMMON_NAME))
myspecies=read.csv(paste(hendl,'dummy.csv',sep=''))%>%
  filter(tolower(Name)%in%Keep.species)
get.this=CAAB%>%
          filter(!SCIENTIFIC_NAME=='Squalus spp.')%>%
          mutate(SCIENTIFIC_NAME=case_when(SCIENTIFIC_NAME=='Squatina australis'~'Squatinidae',  #reset for Families
                                           SCIENTIFIC_NAME=='Pristiophorus cirratus'~'Pristiophoridae',
                                           SCIENTIFIC_NAME=='Squalus megalops'~'Squalus spp.',
                                           SCIENTIFIC_NAME=='Orectolobus maculatus'~'Orectolobidae',
                                           TRUE~SCIENTIFIC_NAME))%>%
          filter(SCIENTIFIC_NAME%in%myspecies$Scientific.name)%>%
            dplyr::select(SPCODE,COMMON_NAME,SCIENTIFIC_NAME)
Store.shp.files=vector('list',length=nrow(get.this))
names(Store.shp.files)=get.this$COMMON_NAME
for(i in 1:length(Store.shp.files))
{
  print(paste('---------i = ',i,'......Downloading shape file for',get.this$SCIENTIFIC_NAME[i]))
  Store.shp.files[[i]]=get_expert_distribution_shp_CAAB(CAAB_species_id=get.this$SPCODE[i],
                                                        spe=get.this$SCIENTIFIC_NAME[i])
}
names(Store.shp.files)=tolower(names(Store.shp.files))
names(Store.shp.files)=case_when(names(Store.shp.files)=='australian angelshark'~'angel sharks',
                                 names(Store.shp.files)=='bronze whaler'~'copper shark',
                                 names(Store.shp.files)=='dusky whaler'~'dusky shark',
                                 names(Store.shp.files)=='greynurse shark'~'grey nurse shark',
                                 names(Store.shp.files)=='spotted wobbegong'~'wobbegongs',
                                 names(Store.shp.files)=='common sawshark'~'sawsharks',
                                 names(Store.shp.files)=='spikey dogfish'~'spurdogs',
                                 TRUE~names(Store.shp.files))
Store.shp.files$hammerheads=Store.shp.files$`smooth hammerhead`


#plot
fn.spatio.temp.catch.dist_old=function(d,Snames,prop.by='fishery',show.prop='blks.with.ktch.over.fished')
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
    mutate(LAT1=abs(LAT1))%>%
    ggplot(aes(LONG1,LAT1))+
    geom_raster(aes(fill = Rel.effort))+
    facet_wrap(~Rango,ncol=2)+
    scale_fill_gradientn(colours=c('ivory2','gold',"red2","darkred"))+
    theme_PA(leg.siz=12.5,axs.t.siz=13,axs.T.siz=19,Sbt.siz=14,strx.siz=14)+
    theme(
      legend.title=element_blank(),
      legend.key.size = unit(.65, 'cm'),
      legend.position = c(0.9, 0.11))+
    ylab(expression('Latitude ('*~degree*S*')'))+xlab(expression('Longitude ('*~degree*E*')'))+
    scale_y_reverse()
  
  
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
    ylab('Proportion of blocks fished with catch')+
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
    scale_y_continuous(sec.axis = sec_axis(~.*coeff, name="Total number of blocks fished"))+
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
fn.spatio.temp.catch.dist=function(d,Snames,FISHRY,core.dist,show.Effort.Catch) 
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
  
  #spatio-temporal effort distribution
  p1=Spatial.effort%>%
    mutate(LAT1=abs(LAT1))%>%
    ggplot(aes(LONG1,LAT1))+
    geom_raster(aes(fill = Rel.effort))+
    facet_wrap(~Rango,ncol=2)+
    scale_fill_gradientn(colours=c('ivory2','gold',"red2","darkred"))+
    theme_PA(leg.siz=12.5,axs.t.siz=13,axs.T.siz=19,Sbt.siz=14,strx.siz=14)+
    theme(
      legend.title=element_blank(),
      legend.key.size = unit(.5, 'cm'),
      legend.position = c(0.92, 0.93))+
    ylab(expression('Latitude ('*~degree*S*')'))+xlab(expression('Longitude ('*~degree*E*')'))+
    scale_y_reverse()
  
  
  #annual blocks with catch over blocks fished within core area
  d11=d1%>%filter(Fishery%in%FISHRY)
  n.spicis=unique(d11$SPECIES)
  d2=vector('list',length(n.spicis))
  names(d2)=n.spicis
  Kr.blks=d2
  for(q in 1:length(n.spicis))
  {
    #1. core area
    core.blocks=d11%>%
      filter(SPECIES==n.spicis[q])%>%
      group_by(BLOCKX)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
      ungroup()%>%
      arrange(-LIVEWT.c)%>%
      mutate(CumKtch=cumsum(LIVEWT.c),
             CumKtch.per=CumKtch/sum(LIVEWT.c,na.rm=T))
    core.blocks=core.blocks[1:which.min(abs(core.blocks$CumKtch.per - core.dist)),]
    core.blocks=core.blocks%>%
      left_join(BLOCKX_lat_long,by='BLOCKX')
    core.blocks=BLOCKX_lat_long%>%
                  filter(LAT<=max(core.blocks$LAT) & LAT>=min(core.blocks$LAT) &
                         LONG<=max(core.blocks$LONG) & LONG>=min(core.blocks$LONG))%>%
                  pull(BLOCKX)
    Kr.blks[[q]]=data.frame(SNAME=capitalize(unique(d11%>%filter(SPECIES==n.spicis[q])%>%pull(SNAME))),
                            BLOCKX=core.blocks)
    
    #2. total number of blocks fished within core area
    Effort_blocks=rbind(Effort.monthly_blocks,Effort.daily_blocks)%>%
      filter(BLOCKX%in%core.blocks)%>%
      count(FINYEAR,BLOCKX)%>%
      group_by(FINYEAR,BLOCKX)%>%
      mutate(n=ifelse(n>0,1,0))%>%
      group_by(FINYEAR)%>%
      summarise(Tot=sum(n))%>%
      data.frame
    Effort.north_blocks=rbind(Effort.monthly.north_blocks,Effort.daily.north_blocks)%>%
      filter(BLOCKX%in%core.blocks)%>%
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
                              arrange(FINYEAR)%>% mutate(Fishery='Northern'))
    Effort_blocks.all=Effort_blocks.all%>%
      left_join(Unik.yr,by='FINYEAR')%>%
      filter(Fishery%in%FISHRY)
    
    #3. positive catch within core area
    N.blocks.with.ktch.yr=d11%>%
      filter(SPECIES==n.spicis[q] & BLOCKX%in%core.blocks)%>%
      count(FINYEAR,SNAME,BLOCKX,Fishery)%>%
      group_by(FINYEAR,SNAME,Fishery)%>%
      mutate(n=ifelse(n>0,1,0))%>%
      group_by(FINYEAR,SNAME,Fishery)%>%
      summarise(n=sum(n,na.rm=T))%>%
      ungroup
    
    #4. Prop of catch out of total blocks fished within core area
    d2[[q]]=left_join(N.blocks.with.ktch.yr,
                      Effort_blocks.all,by=c('FINYEAR','Fishery'))%>%
      mutate(prop=n/Tot)%>%
      group_by(SNAME,Fishery)%>%
      ungroup()                        
  }
  Nfact=2
  Siz=4
  if(length(n.spicis)<5)
  {
    Nfact=1
    Siz=4.5
  }
  d2=do.call(rbind,d2)
  Kr.blks=do.call(rbind,Kr.blks)
  coeff=max(d2$Tot)
  ylim.prim <- c(0, 1.1)   
  ylim.sec <- c(0, max(d2$Tot)) 
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- b*(ylim.prim[1] - ylim.sec[1])

  Sp.order=capitalize(d2%>%group_by(SNAME)%>%summarise(n=sum(n))%>%arrange(-n)%>%pull(SNAME))
  base = d2%>%
    mutate(yr=as.numeric(substr(FINYEAR,1,4)),
           SNAME=capitalize(SNAME),
           SNAME=factor(SNAME,levels=Sp.order))%>%
    ggplot(aes(yr,prop,colour='Proportion with catch'))+
    geom_point(size=2)+
    geom_line()+
    ylab('Proportion of blocks fished with catch within core area')+
    xlab('Financial year')+
    theme_PA(axs.t.siz=13,axs.T.siz=19,Sbt.siz=14,leg.siz=14,strx.siz=16)+
    theme(legend.position = c(0.8, 0.05),  
          legend.title=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=2)))+
    geom_line(aes(x=yr, y = a + Tot*b,color='Number of blocks fished'),size=1.25,alpha=0.35)+
    scale_y_continuous(limits=ylim.prim,
                       sec.axis = sec_axis(~ (. - a)/b, name="Total number of blocks fished within core area"))+
    theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
    theme(axis.title.y.left = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
  main <- base + facet_wrap(~ SNAME,ncol=Nfact)
  
  add.trend=FALSE
  if(add.trend)
  {
    base=base+
      stat_smooth(method = "lm",formula = y ~ x, geom = "smooth",se = FALSE)+
      stat_regline_equation(aes(label=paste(..eq.label..,..rr.label..,sep="~~~~"), color=Fishery),
                            label.y = 1,label.x = quantile(d2$year,.01), size = Siz)
  }
  main <- base + 
    facet_wrap(~ SNAME,ncol=Nfact) +
    scale_color_manual(name = "",
                       breaks = c("Proportion with catch","Number of blocks fished"),
                       values = c("Proportion with catch" = "#F8766D","Number of blocks fished"='black')) 
  
  #dis.neims=sort(unique(d2$SNAME))
  dis.neims=tolower(Sp.order)
  nmax_rep <- length(dis.neims)
  SSIZ=4.5
  vp.wi=0.3
  vp.he = 0.7
  if(nmax_rep>6)
  {
    SSIZ=3.5
    vp.wi=0.4
    vp.he=0.8
  }
  Kr.blks.sorted=Kr.blks%>%mutate(SNAME=factor(SNAME,levels=Sp.order))
  this.shp.file=Store.shp.files[match(dis.neims,names(Store.shp.files))]
  insets <- lapply(seq_len(nmax_rep), function(i) {
    ggplot(ozmap_states) +
      xlim(113,133)+ ylim(-37,-12)+
      geom_sf(fill='orange',alpha=.4)+
      geom_sf(data=sf::st_as_sf(this.shp.file[[match(dis.neims[i],names(this.shp.file))]]),alpha=.4, fill='chartreuse3')+
      geom_point(data=Kr.blks.sorted%>%left_join(BLOCKX_lat_long,by='BLOCKX'),aes(LONG,LAT),alpha=.6,size=.5)+
      theme_bw(base_size=9) +  
      ggforce::facet_wrap_paginate(~ SNAME, nrow = 1, ncol = 1, page = i)+
      guides(fill="none",colour = "none", x = "none", y = "none") +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            axis.title = element_blank(),
            panel.background = element_rect(fill = "transparent",
                                           colour = NA_character_),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent",
                                           colour = NA_character_),
            legend.box.background = element_rect(fill = "transparent"),
            legend.key = element_rect(fill = "transparent"),
            axis.text = element_text(size = 8),
            axis.text.x = element_text(angle = 90),
            axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      annotate(geom="text", x=120, y=-26, label="Core area",color="transparent",size=SSIZ)
  })
  insets2 <- tibble(x = rep(0.01, nmax_rep),
                   y = rep(0.01, nmax_rep),
                   plot = insets,
                   SNAME = factor(Sp.order),levels=Sp.order)
  p2=main +
    geom_plot_npc(data = insets2, 
                  aes(npcx = x, npcy = y, label = plot,
                      vp.width = vp.wi, vp.height = vp.he))
  
  
  if(show.Effort.Catch=='combined')
  {
    figure=ggarrange(plotlist=list(p1,p2),ncol=2,heights=c(1,1))+
      theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
    print(figure)
    ggsave(paste(Rar.path,'/Spatio.temporal.catch_combined_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = wiz,height = heiz,compression = "lzw")
  }
  
  if(show.Effort.Catch=='separated')
  {
    print(p1)
    ggsave(paste(Rar.path,'/Spatio.temporal.effort.tiff',sep=''),width = 6,height = 10,compression = "lzw")
    
    print(p2)
    Wi=7
    Hei=12
    if(length(Snames)<5)
    {
      Wi=7
      Hei=10 
    }
    ggsave(paste(Rar.path,'/Spatio.temporal.catch_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = Wi,height = Hei,compression = "lzw")
    
  }
  
  return(list(Fished.blks.Fishery=Effort_blocks.all,
              prop.ktch_over.fished.bocks=d2))
  rm(Effort_blocks,Effort.north_blocks,NSF.no.fishing,Effort_blocks.all,d2)
}
Store.spatial.temporal.ktch=Lista.sp.outputs
if('additional.sp'%in%names(Store.spatial.temporal.ktch)) Store.spatial.temporal.ktch=Store.spatial.temporal.ktch[-match('additional.sp',names(Store.spatial.temporal.ktch))]
for(l in 1:length(Store.spatial.temporal.ktch))
{
  print(paste("Spatio temporal distribution plots for  ------",names(Lista.sp.outputs)[l]))
  wiz = 12
  heiz = 12
  get.dis.sp=Lista.sp.outputs[[l]]
  if(names(Lista.sp.outputs)[l]=="Other.sp")
  {
    get.dis.sp=c(get.dis.sp,'hammerheads')
    wiz = 13.5
    heiz = 14
  }
  Store.spatial.temporal.ktch[[l]]=fn.spatio.temp.catch.dist(d=Spatio.temp.dat,
                                                             Snames=get.dis.sp,
                                                             FISHRY=c('Southern'),  #c('Northern','Southern')
                                                             core.dist=0.95,
                                                             show.Effort.Catch='separated')   #'combined'
}
clear.log('fn.spatio.temp.catch.dist')
clear.log('Spatio.temp.dat')
clear.log('Effort.daily_blocks')
clear.log('Effort.monthly_blocks')
clear.log('Effort.daily_blocks')
clear.log('Effort.monthly.north_blocks')
clear.log('Effort.daily.north_blocks')
