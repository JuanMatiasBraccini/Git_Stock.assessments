library(dplyr)
Rec.hndl="C:/Matias/Data/Catch and Effort/Recreational/I.Survey."
#Rec.fish.catch.2011.12=read.csv(paste(Rec.hndl,"2011_12.csv",sep=''),stringsAsFactors=F) #Ryan et al 2013 
#Rec.fish.catch.2013.14=read.csv(paste(Rec.hndl,"2013_14.csv",sep=''),stringsAsFactors=F) #Ryan et al 2015 
#Rec.fish.catch.2015.16=read.csv(paste(Rec.hndl,"2015_16.csv",sep=''),stringsAsFactors=F) #Ryan et al 2015 
#Rec.fish.catch=rbind(Rec.fish.catch.2011.12,Rec.fish.catch.2013.14,Rec.fish.catch.2015.16)
Rec.fish.catch=read.csv(paste(Rec.hndl,"csv",sep=''),stringsAsFactors=F)
Rec.fish.catch=Rec.fish.catch%>%mutate(Bioregion=ifelse(Bioregion%in%c("Gasconye","Gascoyne Coast"),"Gascoyne",Bioregion),
                                       FINYEAR=Survey,
                                       Kept.Number=Estimated.Kept.Catch..by.numbers.,
                                       Rel.Number=Estimated.Released.Catch..by.numbers.,
                                       Common.Name=Species.Common.Name,
                                       RSE=Total.RSE....)%>%
                               filter(RSE<=50)%>%  #Karina Ryan recommendation (too inaccurate if higher)  
                        dplyr::select(Common.Name,Kept.Number,Rel.Number,Bioregion,FINYEAR)

AVG.WT=data.frame(Common.Name=c("Bronze Whaler","Greynurse Shark","Gummy Sharks","Scalloped hammerhead",
                                "Smooth hammerhead","Sandbar Shark","Tiger Shark","Whaler Sharks","Wobbegong",             
                                "Other Shark","Western Shovelnose Ray","Rays & Skates","Sawshark",             
                                "Port Jackson Shark","Whiskery Shark","School Shark","Blacktip Reef Shark",   
                                "Dusky Whaler","Lemon Shark","Whitetip Reef Shark"),
                  AVG.wt=c(10,15,5,15,15,10,30,10,10,5,5,5,5,5,5,5,5,15,15,5), #assumed 
                  PCM.rec=c(.1,.1,.2,.4,.4,.1,.1,.1,.1,.1,.1,.1,.1,.05,.2,.2,.1,.1,.1,.1))   #assumed post capture mortality of released sharks        
AVG.WT$Common.Name=as.character(AVG.WT$Common.Name)

#WA population for rec catch recons
WA.population=read.csv("C:/Matias/Data/AusBureauStatistics.csv",stringsAsFactors=F)


#fix species names
Rec.fish.catch=Rec.fish.catch%>%
        mutate(Common.Name=ifelse(Bioregion%in%c("Gascoyne","North Coast") & 
                                         Common.Name=="Bronze Whaler","Dusky Whaler",
                           ifelse(Common.Name%in%c('Whaler Shark – other',"Whaler & Weasel Sharks"),'Whaler Sharks',
                           ifelse(Common.Name=='Hammerhead Shark','Hammerhead Sharks',
                           ifelse(Common.Name%in%c('Other Rays Skates','Unspecified Rays Skates'),'Rays & Skates',
                           ifelse(Common.Name%in%c('Other Sharks','Unspecified Shark',"Sharks"),'Other Shark',
                                  Common.Name))))))

#Reapportion 'whaler sharks' among all reported whaler species
Whaler.prop=Rec.fish.catch%>%
                filter(Common.Name%in%c("Dusky Whaler","Sandbar Shark","Bronze Whaler","Tiger Shark",
                                        "Blacktip Reef Shark","Lemon Shark","Whitetip Reef Shark"))%>%
                group_by(Common.Name,Bioregion,FINYEAR)%>%
                summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
                data.frame
Whaler.prop$Kept.Number.prop=Whaler.prop$Kept.Number/sum(Whaler.prop$Kept.Number)
Whaler.prop$Rel.Number.prop=Whaler.prop$Rel.Number/sum(Whaler.prop$Rel.Number)
Whaler.prop=Whaler.prop%>%dplyr::select(-c(Kept.Number,Rel.Number))

Whaler.Sharks=Rec.fish.catch%>%
                filter(Common.Name=="Whaler Sharks")%>%
                summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
                data.frame
Whaler.reap=cbind(Whaler.prop,Whaler.Sharks)%>%
                mutate(Kept.Number=Kept.Number*Kept.Number.prop,
                       Rel.Number=Rel.Number*Rel.Number.prop)%>%
        dplyr::select(names(Rec.fish.catch))
     
Rec.fish.catch=Rec.fish.catch%>%
                filter(!Common.Name=="Whaler Sharks")
Rec.fish.catch=rbind(Rec.fish.catch,Whaler.reap)

#Reapportion 'other sharks' among all reported shark species                            
Shark.prop=Rec.fish.catch%>%
        filter(!Common.Name%in%c("Other Shark","Western Shovelnose Ray","Rays & Skates","Other Rays and Skates"))%>%
        group_by(Common.Name,Bioregion,FINYEAR)%>%
        summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
        data.frame
Shark.prop$Kept.Number.prop=Shark.prop$Kept.Number/sum(Shark.prop$Kept.Number)
Shark.prop$Rel.Number.prop=Shark.prop$Rel.Number/sum(Shark.prop$Rel.Number)
Shark.prop=Shark.prop%>%dplyr::select(-c(Kept.Number,Rel.Number))


Other.Sharks=Rec.fish.catch%>%
               filter(Common.Name%in%c("Other Shark"))%>%
        summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
        data.frame

Other.reap=cbind(Shark.prop,Other.Sharks)%>%
        mutate(Kept.Number=Kept.Number*Kept.Number.prop,
               Rel.Number=Rel.Number*Rel.Number.prop)%>%
        dplyr::select(names(Rec.fish.catch))

Rec.fish.catch=Rec.fish.catch%>%
        filter(!Common.Name=="Other Shark")
Rec.fish.catch=rbind(Rec.fish.catch,Other.reap)


#Reapportion 'hammerheads'
Rec.fish.catch=Rec.fish.catch%>%
                mutate(Common.Name=ifelse(Common.Name=='Hammerhead Sharks' &
                                          Bioregion%in%c("Gascoyne","North Coast"),"Scalloped hammerhead",
                                   ifelse(Common.Name=='Hammerhead Sharks' &
                                          Bioregion%in%c("West Coast","South Coast"),"Smooth hammerhead",
                                   Common.Name)))

fn.rec=function(DAT)
{
        AGG=DAT%>%left_join(AVG.WT,by="Common.Name")%>%
                  mutate(LIVEWT.c=round(Kept.Number+Rel.Number*PCM.rec)*AVG.wt)%>%
                group_by(FINYEAR,Common.Name,Bioregion)%>%
                summarise(LIVEWT.c=sum(LIVEWT.c))%>%
                data.frame
        LisT=vector('list',length(unique(AGG$Common.Name)))
        names(LisT)=unique(AGG$Common.Name)
        for(i in 1:length(LisT))
        {
                LisT[[i]]=subset(AGG,Common.Name==names(LisT)[i],select=c(Common.Name,LIVEWT.c,FINYEAR,Bioregion))
        }
        return(LisT)
}
Rec.ktch=fn.rec(Rec.fish.catch)  
Part.rate.89=26.6  #Ryan et al 2012
Part.rate.00=28.5
Part.rate=lm(c(Part.rate.89,Part.rate.00)~c(1989,2000))
Part.rate.pred=coef(Part.rate)[1]+coef(Part.rate)[2]*WA.population$Year
Part.rate.pred=(Part.rate.pred/100)*WA.population$Population
Part.rate.pred=Part.rate.pred/Part.rate.pred[match(2011,WA.population$Year)]
Part.rate.pred=data.frame(P.rate=Part.rate.pred[1:(length(Part.rate.pred)-1)],
                          FinYear=paste(WA.population$Year[1:(length(WA.population$Year)-1)],"-",
                                        substr(WA.population$Year[2:length(WA.population$Year)],start=3,stop=4),sep=""))

back.fill=function(dat)
{
        id=match(unique(dat$FinYear),Part.rate.pred$FinYear)
        Regns=unique(dat$Bioregion)
        Dummy=vector('list',length(Regns))
        for (d in 1:length(Dummy))
        {
                a=subset(dat,Bioregion==Regns[d])
                Mat=matrix(NA,nrow=length(Part.rate.pred$FinYear),ncol=2)
                Mat[,2]=mean(a$LIVEWT.c)*Part.rate.pred$P.rate
                Mat=as.data.frame(Mat)
                names(Mat)=c("FINYEAR","LIVEWT.c")
                Mat$FINYEAR=Part.rate.pred$FinYear
                Mat$zone=Regns[d]
                Mat$Common.Name=dat$Common.Name[1]
                Dummy[[d]]=Mat
        }
        Dummy=do.call(rbind,Dummy)
        return(Dummy)
}
for(i in 1:length(Rec.ktch)) Rec.ktch[[i]]=back.fill(Rec.ktch[[i]])
Rec.ktch=do.call(rbind,Rec.ktch)
row.names(Rec.ktch)=NULL

#Change Backtip reef shark to spinner if in the south
Rec.ktch=Rec.ktch%>%mutate(Common.Name=ifelse(Common.Name=="Blacktip Reef Shark" & 
                        zone%in%c('South Coast','West Coast'),'Spinner Shark',Common.Name))