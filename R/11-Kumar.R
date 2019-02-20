#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 11-Kumar 1997
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

###Inputs:
#  Esc=1353   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  press      [mb]        (local barometric)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceKumar<-function(){
	#extraterrestrial irradiance
	Esc=1353
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
	#absolute air mass
	mk=((1229+(614*cos(sza))^2)^0.5-614*cos(sza))*press/1013.25
	
	#direct beam irradiance
	Ebnkumar=Eext*0.56*(exp(-0.65*mk)+exp(-0.095*mk))
	#diffuse horizontal irradiance
	Edhkumar=(0.2710*Eext-0.294*Ebnkumar)*cos(sza)
	#global horizontal irradiance
    Eghkumar=Ebnkumar*cos(sza)+Edhkumar
	
	 ###Quality control
     lower=0
     Ebnkumar[Ebnkumar<lower]=0
     Edhkumar[Edhkumar<lower]=0
     Eghkumar[Eghkumar<lower]=0
	 return(list(Ebnkumar, Edhkumar, Eghkumar))
     }


