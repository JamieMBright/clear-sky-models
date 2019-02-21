###A simple example

###Inputs data example is in 'adelaide_airport_inputs_sample.csv', it's one day (20150119 ) data of adelaide airport

###Inputs LIST from 'adelaide_airport_inputs_sample.csv' (some may not be used):
###Year,Month,Day,Dayth,Hour,Minute,Second,ghi,dni,dif,sza,NO2,press,albedo,ang_alpha,ang_beta,aod550,aod700,broadbandaod,ozone,wv,Tem,Tdpt,aerosol_scattering,TL2R,TL2D,TL2I,TL2Gu,TL2M,TL2Gr

###please copy the following line to the R console

#         source(file='Example.R')

###please copy the upper line to the R console
###make sure 'Example.R','xx-model.R'  and 'adelaide_airport_inputs_sample.csv' are in the working folder

###Codes
#load all inputs data
adelaide_airport=read.table('adelaide_airport_inputs_sample.csv',header=TRUE,sep=',')
attach(adelaide_airport)

#load clear sky model function
source(file='6-ASHRAE.R')
#you can change 'xx-model.R' to test other models

#Run function
#Because model TJ can produce 3 components Ebn,Edh,Egh, the function will return a list with 3 components 
Ebn=IrradianceASHRAE()[[1]]
Edh=IrradianceASHRAE()[[2]]
Egh=IrradianceASHRAE()[[3]]
#you can change IrradianceModel() from 'xx-model.R' to test other models

#Visualization
plot(Ebn,type="l",lty=2,ylim=c(0,1200))
lines(Edh,col="red",lty=3)
lines(Egh,col="blue",lty=4)
lines(ghi,col="yellow",lty=5)
title("comparison",cex.main=.8)
legend("topleft",inset=.05,c("Ebn","Edh","Egh","ghi"), cex=.5,lty=c(2,3,4,5),col=c("black","red","blue","yellow"))