#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 14-Atwater and Ball_1 1981
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Atwater, M. A., & Ball, J. T. (1981). A surface solar radiation model for cloudy atmospheres. Monthly weather review, 109(4), 878-888.

###Inputs:
#  Esc=1353   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  press      [mb]        (local barometric)
#  wv         [atm.cm]    (total columnar amount of water vapour)

###Outputs:
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceAtwater1<-function(){	
	#extraterrestrial irradiance
	Esc=1353
    Eext=Esc*(1+0.033*cos(2*pi*Dayth/365))
    
	# Air Mass for Atwater	
    am=35/(1224*(cos(sza))^2+1)^0.5 
    ama=am*press/1013.25
    
    #the transmittance after Rayleigh scattering TR and absorptance of permanent gases Tg 
    presskpa=press/10    #our pressure unit is hpa/mb, in paper the unit is kpa. so a unit conversion is applied
    TRTg=1.021-0.0824*((949*presskpa*(10^-5)+0.051)*am)^0.5
    
    #water vapor transmittance
    uW=wv*am
    TW=1-2.9*uW/(5.925*uW+(1+141.5*uW)^0.635)
    
    #the transmittance function for multiple layers of clouds
    Tc=1 #we calculate clear sky irraiance, so Tc=1
    
    #the transmittance after absorptance and reflectance of  aerosols
    Tp=1 #author assumes Tp=1
    
    #atmospheric reflectance
    rs=0.0685  #as suggested in ATWATER1980, originally by Lacis and Hansen <A Parameterization for Absorption of Solar Radiation in the Earth's Atmosphere>
    #ground albedo
    rg=-0.0139+0.0467*tan(sza)      #author limits 1>=rg>=0.03
    rg[which(rg>1)]=1
    rg[which(rg<0.03)]=0.03
    
    #global horizontal irradiance
    Eghatwater1 = Eext*cos(sza)*TRTg*TW*Tc*Tp/(1-rg*rs)
    
    ###Quality control
    lower=0
    Eghatwater1[Eghatwater1<lower]=0
    return(Eghatwater1)
    }
    
