#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 70-Bird 1981
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Bird, R. E., & Hulstrom, R. L. (1981). Simplified clear sky model for direct and diffuse insolation on horizontal surfaces (No. SERI/TR-642-761). Solar Energy Research Inst., Golden, CO (USA).

###Inputs:
#  Esc=1353     [Wm-2]             (Solar constant)
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
IrradianceBird<-function(){
    #extraterrestrial irradiance
	Esc=1353
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eex=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
	
	#Air Mass 
	am=(cos(sza)+0.15*((pi/2-sza)/pi*180+3.885)^-1.25)^-1
    ama=am*press/1013.25
    
    #transmittaance of Rayleigh scattering
    TR=exp(-0.0903*(ama^0.84)*(1+ama-ama^1.01))
    
    #transmittaance of Ozone absorptance
    XO=ozone*am
    TO=1-0.1611*XO*(1+139.48*XO)^(-0.3035)-0.002715*XO*(1+0.044*XO+0.0003*XO^2)^-1	
    
    #transmittaance of absorptance by uniformly mixed gases(CO2,O2) 
    TUM=exp(-0.0127*ama^0.26)
    
    #transmittaance of water vapor absorptance
    XW=wv*am
    TW=1-2.4959*XW*((1+79.034*XW)^0.6828+6.385*XW)^-1

    #transmittance of aerosol absorptance and scattering
    broadaod=ang_beta*(0.2758*0.38^-ang_alpha+0.35*0.5^-ang_alpha)
    TA=exp((-broadaod^0.873)*(1+broadaod-broadaod^0.7088)*am^0.9108) 

    #direct beam irradiance
    Ebnbird=Eex*0.9662*TA*TW*TUM*TO*TR
	
    #transmittance of aerosol absorptance
    K1=0.1 #as suggested by author
    TAA=1-K1*(1-am+am^1.06)*(1-TA)

    #transmittance of aerosol scattering
    TAS=TA/TAA
     
    BA=0.84 #as suggested by author
    Eas=Eex*cos(sza)*0.79*TO*TUM*TW*TAA*(0.5*(1-TR)+BA*(1-TAS))/(1-am+am^1.02)
    
    rs=0.0685+(1-BA)*(1-TAS)  #sky albedo
    rs[which(rs<0)]=0
    rs[which(rs>1)]=1
    rg=albedo  #ground albedo
    
    #global horizontal irradiance
    Eghbird=(Ebnbird*cos(sza)+Eas)/(1-rg*rs)
    #diffuse horizontal irradiance
    Edhbird=Eghbird-Ebnbird*cos(sza)

	###Quality control
    lower=0
    Ebnbird[Ebnbird<lower]=0
    Edhbird[Edhbird<lower]=0
    Eghbird[Eghbird<lower]=0
	return(list(Ebnbird, Edhbird, Eghbird))
	}



