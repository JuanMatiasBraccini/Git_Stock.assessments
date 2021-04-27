library(tidyverse)
library(ggplot2)

# DATA SECTION ------------------------------------------------------------
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')
All.species.names=read.csv(handl_OneDrive("Data/Species_names_shark.only.csv"),stringsAsFactors = F)

fn.in=function(NM) read.csv(paste(handl_OneDrive('Analyses/Data_outs/'),NM,sep=""),stringsAsFactors = F)


#2.1 Catch_WA Fisheries

  #Historic
Hist.expnd=fn.in(NM='recons_Hist.expnd.csv')

  #Ammended reported catch including discards
Data.monthly=fn.in(NM='recons_Data.monthly.csv')
Data.monthly.north=fn.in(NM='recons_Data.monthly.north.csv')

  #TEPS
Greynurse.ktch=fn.in(NM='recons_Greynurse.ktch.csv')
TEPS_dusky=fn.in(NM='recons_TEPS_dusky.csv')

  #Droplines Western Rock Lobster
WRL.ktch=fn.in(NM='Wetline_rocklobster.csv')


#2.2. Catch of non WA Fisheries

  #Taiwanese gillnet and longline
Taiwan.gillnet.ktch=fn.in(NM='recons_Taiwan.gillnet.ktch.csv')
Taiwan.longline.ktch=fn.in(NM='recons_Taiwan.longline.ktch.csv')

  #Indonesian illegal fishing in Australia waters
Indo_total.annual.ktch=fn.in(NM='recons_Indo.IUU.csv') 

  #AFMA's GAB & SBT fisheries
GAB.trawl_catch=fn.in(NM='recons_GAB.trawl_catch.csv') 
WTBF_catch=fn.in(NM='recons_WTBF_catch.csv') 

  #SA Marine Scalefish fishery
Whaler_SA=fn.in(NM='recons_Whaler_SA.csv') 



#2. Reconstructed recreational catch time series
Rec.ktch=fn.in(NM='recons_recreational.csv')


# PROCEDURE SECTION ------------------------------------------------------------

All.species.names=All.species.names%>%
  mutate(Name=tolower(Name))%>%
  rename(SNAME=Name)

# Manipulate commercial catch
fn.sel.columns=function(d) d%>%dplyr::select(names(Hist.expnd))

Data.monthly=fn.sel.columns(Data.monthly)
Data.monthly.north=fn.sel.columns(Data.monthly.north)
Greynurse.ktch=fn.sel.columns(Greynurse.ktch)
TEPS_dusky=fn.sel.columns(TEPS_dusky)
WRL.ktch=fn.sel.columns(WRL.ktch)
Taiwan.gillnet.ktch=fn.sel.columns(Taiwan.gillnet.ktch)
Taiwan.longline.ktch=fn.sel.columns(Taiwan.longline.ktch)
Indo_total.annual.ktch=fn.sel.columns(Indo_total.annual.ktch)
GAB.trawl_catch=fn.sel.columns(GAB.trawl_catch)
WTBF_catch=fn.sel.columns(WTBF_catch)
Whaler_SA=fn.sel.columns(Whaler_SA)

Com.ktch=rbind(Hist.expnd,
               Data.monthly,
               Data.monthly.north,
               Greynurse.ktch,
               TEPS_dusky,
               WRL.ktch,
               Taiwan.gillnet.ktch,
               Taiwan.longline.ktch,
               Indo_total.annual.ktch,
               GAB.trawl_catch,
               WTBF_catch,
               Whaler_SA)%>%
          group_by(FINYEAR,SPECIES)%>%
          summarise(LIVEWT.c=sum(LIVEWT.c))%>%
          left_join(All.species.names,by='SPECIES')


#Manipulate recreational catch
Rec.ktch=Rec.ktch%>%
            group_by(FINYEAR,Common.Name)%>%
            summarise(LIVEWT.c=sum(LIVEWT.c))%>%
            mutate(Common.Name=ifelse(Common.Name=="dogfishes","spurdogs",
                               ifelse(Common.Name=="greynurse shark","grey nurse shark",
                               ifelse(Common.Name=="thresher shark","thresher sharks",
                               ifelse(Common.Name=="bronze whaler","copper shark",
                               Common.Name)))))%>%
            left_join(All.species.names,by=c('Common.Name'='SNAME'))%>%
            mutate(SNAME=tolower(Common.Name))%>%
            dplyr::select(names(Com.ktch))


#Combine commercial and recreational
Tot.ktch=rbind(Com.ktch,Rec.ktch)%>%
            group_by(FINYEAR,SNAME)%>% 
            summarise(LIVEWT.c=sum(LIVEWT.c)/1000)%>%    #in tonnes
            rename(LIVEWT_tonnes=LIVEWT.c)%>%
            mutate(SNAME=ifelse(SNAME=="angel shark","angel sharks",
                         ifelse(SNAME=="dusky whaler","dusky shark",
                         ifelse(SNAME=="grey nurse shark","greynurse shark",
                         ifelse(SNAME=="gummy sharks","gummy shark",
                         ifelse(SNAME=="tawny shark","tawny nurse shark",
                         ifelse(SNAME=="thresher shark","thresher sharks",
                         ifelse(SNAME=="spotted wobbegong","wobbegongs",
                         ifelse(SNAME%in%c("southern sawshark","common sawshark"),"sawsharks",
                                SNAME)))))))))

               

#Data for Northern Territory assessment of spot tail, australian blacktip, common blacktip sharks
NT.blktips=c("spot-tail shark","blacktips")
For.NT_blktips=Tot.ktch%>%filter(SNAME%in%NT.blktips)
#ggplot(For.NT_blktips)+geom_line(aes(x=as.numeric(substr(FINYEAR,1,4)),y=LIVEWT_tonnes,colour =SNAME))
write.csv(For.NT_blktips,handl_OneDrive("Analyses/Catch and effort/Data_Resquests/Northern_Territory/Blacktips.and.spot.tail.total.reconstructed.catches.csv"),row.names = F)
