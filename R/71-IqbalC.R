#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 71-Iqbal-C 1983
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Iqbal, M. (2012). An introduction to solar radiation. Elsevier.

###Inputs:
#  Esc=1367     [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  press        [mb]               (local barometric)
#  albedo       [double]           (surface/ground albedo)
#  ang_alpha    [dimensionless]    (Angstrom_exponent, also known as alpha)
#  ang_beta     [dimensionless]    (Angstrom turbidity coefficient)
#  ozone        [atm.cm]           (total columnar amount of ozone)
#  wv           [atm.cm]           (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceIqbalc<-function(){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))

    #air mass
    am=(cos(sza)+0.15*(3.885+(pi/2-sza)/pi*180)^-1.253)^-1
    ama=am*press/1013.25
   
    #the transmittance by Rayleigh scattering
    TR=exp(-0.0903*(ama)^0.84*(1+ama-ama^1.01)) 

     #the transmittance by ozone
     U3=ozone*am       
     alpho=(0.1611*U3*(1+139.48*U3)^-0.3035 - 0.002715*U3*(1.0+0.044*U3+0.0003*U3^2)^-1)
     TO=1-alpho

     #the transmittance by uniformly mixed gases
     TG=exp(-0.0127*ama^0.26)

     #the transmittance by water vapor
     PW=am*wv*(press/1013.25)^0.75
     alphw=2.4959*PW*((1+79.034*PW)^0.6828+6.385*PW)^-1
     TW=1-alphw

     #the transmittance by Aerosol
     kaa=ang_beta*(0.2758*0.38^-ang_alpha+0.35*0.5^-ang_alpha)     
     TA=exp(-kaa^0.873*(1+kaa-kaa^0.7088)*ama^0.9108)

     #direct beam irradiance
     Ebniqbalc=0.9751*Eext*TR*TO*TG*TW*TA

     #the transmittance of direct radiation due to aerosol absorptance
     w0=0.9  #as suggested by author from Bird Hulstrom
     TAA=1-(1-w0)*(1-ama+ama^1.06)*(1-TA)
     #the Rayleigh-scattered diffuse irradiance after the first pass through the atmosphere
     Edr = 0.79*Eext*cos(sza)*TO*TG*TW*TAA*0.5*(1-TR)/(1-ama+ama^1.02)
     
     #the aerosol scattered diffuse irradiance after the first pass through the atmosphere
     TAS=TA/TAA
     Fc=0.84   #as suggested by author
     Eda=0.79*Eext*cos(sza)*TO*TG*TW*TAA*Fc*(1-TAS)/(1-ama+ama^1.02)

     #global horizontal irradiance
     rg=albedo   #ground albedo
     ra=0.0685+(1-Fc)*(1-TAS)
     Eghiqbalc=(Ebniqbalc*cos(sza)+Edr+Eda)/(1-rg*ra)
     #diffuse horizontal irradiance
     Edhiqbalc=Eghiqbalc-Ebniqbalc*cos(sza)
     
     ###Quality control
     lower=0
     Ebniqbalc[Ebniqbalc<lower]=0
     Edhiqbalc[Edhiqbalc<lower]=0
     Eghiqbalc[Eghiqbalc<lower]=0
	 return(list(Ebniqbalc, Edhiqbalc, Eghiqbalc))
	 }


