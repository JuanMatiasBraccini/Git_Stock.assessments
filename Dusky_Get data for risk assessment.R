#script for creating input data for Risk Assessments

#Steps: Section A: Bring in data (this brings in all the data series from "Organise data.R")

rm(list=ls(all=TRUE))


# ................. Section A: BRING IN INPUT DATA ......................
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Population dynamics/Organise data.R")

#Select species and year of assessment
SPEC=18003     #note: Dusky includes C. brachyurus as well (in "Organise data.R")
Spec="Dusky"
species="BW"
AssessYr=2017
TL.bins.cm=5
Data.yr="2015-16"         #last year of catch
Frst.yr.ktch="1975-76"    #first year of catch
#Minimun number of samples of size composition
Min.obs=10
Min.shts=10


#select whether to use all catch series or only WA
#Ktch.source="ALL"
Ktch.source="WA.only"

#Create data inputs 
  #define is using conventional tagging and effort in model
add.conv.tag="YES"
Add.Effort="NO"
fn.input.data(SP=species,Yr.assess=AssessYr,Conv.cal.mn.to.fin.mn="NO",           
              Historic.Ktch="NO",Bin.size=TL.bins.cm,What.Efrt="km.gn.days")