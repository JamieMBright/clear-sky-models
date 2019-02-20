#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 07-Sharma 1965
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Sharma, M. R., & Pal, R. S. (1965). Interrelationships between total, direct, and diffuse solar radiation in the tropics. Solar Energy, 9(4), 183-192.
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  Esc=1366.1   [Wm-2]     (Solar constant)
#  sza          [radians]  (zenith_angle) 

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceSharma<-function(){
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
	B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
	#direct normal irradiance
	EbnSharma=1.842*(Eext/2)*cos(sza)/(0.3135+cos(sza)) 
	#global horizontal irradiance
	EghSharma=4.5*(Eext/(2*60))+1.071*EbnSharma*cos(sza) 
	#diffuse horizontal irradiance
	EdhSharma=EghSharma-EbnSharma*cos(sza)
	
    ###Quality control
    lower=0
    EbnSharma[EbnSharma<lower]=0
    EdhSharma[EdhSharma<lower]=0
    EghSharma[EghSharma<lower]=0
	return(list(EbnSharma, EdhSharma, EghSharma))
	}


