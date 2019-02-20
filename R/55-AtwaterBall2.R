#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 55-Atwater and Ball_2 1981 
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Bird, R., & Hulstrom, R. (1981). A simplified clear sky model for direct and diffuse insolation on horizontal surfaces. SERI. Solar Energy Research Institute, Golden, Colorado, 642-661.

###Inputs:
#  Esc=1353         [Wm-2]              (Solar constant)
#  sza              [radians]           (zenith_angle) 
#  press            [mb]                (local barometric)
#  albedo           [double]            (surface/ground albedo)
#  Broadband_AOD    [dimensionless]     (broadband aerosol optical depth)
#  wv               [atm.cm]            (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceAtwater2<-function(){	
	#extraterrestrial irradiance
	Esc=1353
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
	
	#Air Mass for Atwater	
    am=35/(1224*(cos(sza))^2+1)^0.5 
    ama=am*press/1013.25
    
    #Rayleigh Transmittance
    pressmilibar=press
    TM=1.021-0.0824*((949*pressmilibar*10^-6+0.051)*am)^0.5
    TMd=1.041-0.16*((949*pressmilibar*10^-6+0.051)*am)^0.5
    
    #Water Vapor Transmittance
    PW=wv
    AW=0.077*(PW*am)^0.3
    
    #Aerosol Optical Depth
    TA=exp(-broadbandaod*ama)
    
    rs=0.0685  #as suggested by author,  originally in Lacis and Hansen "A Parameterization for Absorption of Solar Radiation in the Earth's Atmosphere "
    rg =albedo   #surface reflectance
     
    #direct normal radiation
	Ebnatwater2=Eext*(TMd-AW)*TA
	
	#global horizontal radiation
    Eghatwater2 = Eext*cos(sza)*(TM-AW)*TA/(1-rg*rs)
    
	#diffuse horizontal radiation
	Edhatwater2=Eghatwater2-Ebnatwater2*cos(sza)    
	
    ###Quality control
    lower=0
    Ebnatwater2[Ebnatwater2<lower]=0
    Edhatwater2[Edhatwater2<lower]=0
    Eghatwater2[Eghatwater2<lower]=0
	return(list(Ebnatwater2,Edhatwater2,Eghatwater2))
	}


