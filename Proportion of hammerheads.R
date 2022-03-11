#Define user
User="Matias"

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')


#Sharks data base 
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))
setwd(handl_OneDrive('Analyses/Hammerhead assessments'))


# Hammerhead species composition by zone and time
HH=DATA%>%
  filter(CAAB_code>19000 & CAAB_code<19010)%>%
  mutate(Lat.band=round(abs(Mid.Lat)))

d1=HH%>%
  filter(!is.na(zone))%>%
  mutate(zone=case_when(zone=="Closed"~ "Ningaloo closure",
                   zone=="Joint" ~"JA NSF",
                   zone=='North' ~ "WA NCSF",
                   TRUE~zone))%>%
  group_by(COMMON_NAME,zone,year)%>%
  summarise(n=sum(Numbers))%>%
  ungroup()%>%
  group_by(zone,year)%>%
  mutate(Proportion = n / sum(n))


#proportions with number of HH individuals sampled
d1%>%
  ggplot(aes(year,Proportion,color=COMMON_NAME))+
  geom_point(aes(size=n))+
  facet_wrap(~zone)+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab('Proportion of hammerhead catch')+
  xlab("Year")


#proportions with number of shots
d1.shots=HH%>%
  distinct(SHEET_NO,.keep_all=T)%>%
  mutate(zone=case_when(zone=="Closed"~ "Ningaloo closure",
                        zone=="Joint" ~"JA NSF",
                        zone=='North' ~ "WA NCSF",
                        TRUE~zone))%>%
  group_by(zone,year)%>%
  summarise(N.shots=n())


d1%>%
  left_join(d1.shots,by=c('zone','year'))%>%
  ggplot(aes(year,Proportion,color=COMMON_NAME))+
  geom_point(aes(size=N.shots))+
  facet_wrap(~zone)+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab('Proportion of hammerhead catch')+
  xlab("Year")
ggsave('Proportion hh with number of shots.tiff', width = 12,height = 8, dpi = 300, compression = "lzw")


d2=HH%>%
  filter(zone%in%c('West','Zone1','Zone2'))%>%
  group_by(COMMON_NAME,year)%>%
  summarise(n=sum(Numbers))%>%
  ungroup()%>%
  group_by(year)%>%
  mutate(Proportion = n / sum(n))

d2%>%
  left_join(HH%>%
              distinct(SHEET_NO,.keep_all=T)%>%
              filter(zone%in%c('West','Zone1','Zone2'))%>%
              group_by(year)%>%
              summarise(N.shots=n())
            ,by=c('year'))%>%
  ggplot(aes(year,Proportion,color=COMMON_NAME))+
  geom_point(aes(size=N.shots))+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab('Proportion of hammerhead catch')+
  xlab("Year")
ggsave('Proportion hh with number of shots combined WC_SC.tiff', width = 12,height = 8, dpi = 300, compression = "lzw")


