#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 13-Fu and Rich 1999 (Release v9.2 2008)
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

###Inputs:
#  Esc=1367   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  altitude   [m]         (location's altitude)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceFurich<-function(altitude){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))  
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
	#the air mass corrected for elevation
	mf=exp(-0.000118*altitude-1.638e-9*altitude^2)/cos(sza)
	#a bulk atmospheric transmittance
	Tb=0.5  #recommended value in Gueymard 2012
	
	#direct beam irradiance
	Ebnfurich=Eext*(Tb^mf)
	
	#diffuse fraction of global normal irradiance 
	P=0.3  #recommended value in Gueymard 2012
	
	#diffuse horizontal irradiance
	Edhfurich=Ebnfurich*(P/(1-P))*cos(sza)
	
	#global horizontal irradiance
	Eghfurich=Ebnfurich*cos(sza)+Edhfurich
	
	###Quality control
	lower=0
    Ebnfurich[Ebnfurich<lower]=0
    Edhfurich[Edhfurich<lower]=0
   	Eghfurich[Eghfurich<lower]=0
    return(list(Ebnfurich, Edhfurich, Eghfurich))
    }


