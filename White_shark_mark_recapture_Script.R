#MArk recapture analysis of WA white sharks

#notes: Three types of tagging types: 1) internal acoustic & plastic fin tag
#                                     2) internal acoustic only (2012, same week)
#                                     3) external acoustic only    


library(dplyr)
setwd('C:\\Matias\\Analyses\\Population dynamics\\White shark\\Mark_recapture')
Data=read.csv('Data.csv',stringsAsFactors = F)


#Recaptures
#Comments from Silas on 25/10/2018:
#     Hamish was originally tagged in SA, we retagged him in Albany.
#     Cupcake was given an external when she swam past the boat one day, we 
#     caught her days later in the same place and gave her an internal.
#     Lyn was tagged in Cockburn Sound in 2012, caught and killed in Esperance 2014

#     Total sharks tagged in WA is 87. 
#     There should be more tag records for both Tania and Bethwyn. 
#     Tania should have records for 2012 and 2013. Bethwyn 2013 and 2014.

#     Any without a release location but with a length and realistic release date 
#     will be from SA. Ours will have the full details. Some of the others that 
#     have no length, sex or release info are either not deployed, lost or we have
#     no idea whether they've been deployed.


#Add 'Lyn'. Tagged and then killed at Wyllies Beach
#note: Dani's email
# "Hi Matias, 
#This is what I found so far: M:\Fisheries Research\FinFish\Shark\Samples to other researchers 
#File name "WP Genetics Info (RussB) Feb 2016" sample numbers R3000/100 and R3000/101 you have the size for the 2 animals captured and killed back in 2014
#Cheers,
#Dani
#"

Killed.dissected="Lyn"

Data.Lyn.recp=Data[match("Lyn",Data$Name2),] #recapture info for Lyn
Data.Lyn.recp[,]=NA
Data.Lyn.recp$Species2="white"
Data.Lyn.recp$Sex2="F"
Data.Lyn.recp$Name2= "Lyn"                       
Data.Lyn.recp$ReleaseDate2="2-Oct-14"     
Data.Lyn.recp$ReleaseLength=3.04
Data.Lyn.recp$ReleaseSite2="Wylie beach"
Data.Lyn.recp$ReleaseLongitude2=121.9971
Data.Lyn.recp$ReleaseLatitude2=33.8325
Data=rbind(Data,Data.Lyn.recp)  


#Add sharks killed as part of imminent threat drum lines    #MISSING: Falcons
Killed=c("Wyllie.killed1","Falcon.killed")
Add.killed=Data[1:length(Killed),]
Add.killed[,]=NA
Add.killed$Species2="white"
Add.killed$Name2=Killed
Add.killed$Sex2=c("F","U")
Add.killed$ReleaseDate2=c("2-Oct-14","1-Jun-16")  
Add.killed$ReleaseLength=c(2.37,4.2)
Add.killed$ReleaseSite2=c("Wylie beach","Falcon Beach")
Add.killed$ReleaseLongitude2=c(121.9971,115.659)
Add.killed$ReleaseLatitude2=c(33.8325,32.5769)
Data=rbind(Data,Add.killed) 


#Known recaptures (Silas pers comm)
Reported.recaptures=c("Tania", "Bethwyn", "Max",Killed.dissected)


Data=Data%>%mutate(
  ReleaseDate=as.Date(ReleaseDate2, "%d-%b-%y"),
  Length=gsub("[^0-9.]", "",  Data$ReleaseLength),
  Latitude=-abs(as.numeric(ReleaseLatitude2)),
  Longitude=as.numeric(ReleaseLongitude2),
  Latitude=ifelse(ReleaseSite2=="Albany",-35.0275,
          ifelse(ReleaseSite2=="Doubtful Islands",-34.3751,
          ifelse(ReleaseSite2=="Alexander Point",-33.8839,
          ifelse(ReleaseSite2=="Salisbury Island",-34.3597,
          Latitude)))),
  Longitude=ifelse(ReleaseSite2=="Albany",117.8840,
          ifelse(ReleaseSite2=="Doubtful Islands",119.6072,
          ifelse(ReleaseSite2=="Alexander Point",122.7633,
          ifelse(ReleaseSite2=="Salisbury Island",123.5532,
          Longitude)))),
  Jurisdiction=ifelse(Longitude<=129,"WA",
          ifelse(Longitude>129 & Longitude<141,"SA",
          ifelse(Longitude>141 & Longitude<150,"Vic",
          ifelse(Longitude>150,"NSW",
          ifelse(grepl('NSW',Name2),"NSW",
          NA)))))
          )

Data$Jurisdiction=with(Data,
    ifelse(is.na(Jurisdiction)& ReleaseSite2%in%
        c("Boulders Beach, Ballina","Evans Head","Hawk's Nest","hawks Nest",
          "Hawks Nest","Lennox head","Lighthouse Beach, Ballina",
          "Lighthouse Head, Ballina","Main Beach, Evans Head",
          "Manning River, Croki","NSW","Rabbit Island","Shelly Beach, Ballina"),"NSW",
    ifelse(is.na(Jurisdiction)& ReleaseSite2%in%
           c("North Neptunes","North Neptune","South Neptune"),"SA",Jurisdiction)))

#Keep only WA tags
WA.tags=subset(Data,Jurisdiction=="WA")


#Identify and remove double tags
WA.tags$double.tag=with(WA.tags,paste(ReleaseDate,Name2,Sex2,Length))
Double.tagged=WA.tags[duplicated(WA.tags$double.tag),]
Double.tagged=subset(WA.tags,double.tag%in%Double.tagged$double.tag)
Double.tagged=Double.tagged[order(Double.tagged$ReleaseDate,Double.tagged$Name2),]
write.csv(Double.tagged,"Double.tagged.csv")
WA.tags=WA.tags[!duplicated(WA.tags$double.tag),]

#Identify recaptures
WA.tags$Recapture=with(WA.tags,ifelse(Name2%in%Reported.recaptures,1,0))



#Aca
#Add type of tag (acoustic and fin tag)
#$TagType

#Select variables of interest
WA.tags=subset(WA.tags,select=c(Code2,Sex2,Name2,ReleaseDate,ReleaseSite2,
                                Length,Latitude,Longitude,Recapture,TagType))