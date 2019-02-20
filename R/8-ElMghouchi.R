#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 08-El Mghouchi 2014
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#El Mghouchi, Y., El Bouardi, A., Choulli, Z., & Ajzoul, T. (2014). New model to estimate and evaluate the solar radiation. International Journal of Sustainable Built Environment, 3(2), 225-234.

###Inputs:
#  Esc=1367   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceElMghouchi<-function(){
	#extraterrestrial irradiance
	Esc=1367
    Eext=Esc*(1+0.034*cos(Dayth-2))
    
    #turbidity atmospheric factor for clear skies 
    TFactor=0.796-0.01*sin(0.986*(Dayth+284))
    
	#direct normal irradiance
	EbnElMghouchi=Eext*TFactor*exp(-0.13/cos(sza))
	
	#diffuse horizontal irradiance
	EdhElMghouchi=120*TFactor*exp(-1/(0.4511+cos(sza)))
	
	#global horizontal irradiance
	EghElMghouchi=EbnElMghouchi*cos(sza)+EdhElMghouchi
	
    ###Quality control
    lower=0
    EbnElMghouchi[EbnElMghouchi<lower]=0
    EdhElMghouchi[EdhElMghouchi<lower]=0
    EghElMghouchi[EghElMghouchi<lower]=0
	return(list(EbnElMghouchi, EdhElMghouchi, EghElMghouchi))
	}
	
	