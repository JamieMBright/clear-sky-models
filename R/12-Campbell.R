#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 12-Campbell 1998
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Campbell, G. S., & Norman, J. M. (1998). An introduction to environmental biophysics. An introduction to environmental biophysics., (Ed. 2).

###Inputs:
#  Esc=1366.1   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  press      [mb]        (local barometric)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceCampbell<-function(){
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
	B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #air mass
	am=press/1013*(cos(sza))^-1
	
	#direct normal irradiance
	tau0=0.7  #Author stated 0.7 to be typical of clear sky conditions. 
	EbnCampbell=Eext*tau0^am
	#diffuse horizontal irradiance
	EdhCampbell=0.3*Eext*cos(sza)*(1-tau0^am)
	#global horizontal irradiance
	EghCampbell=EbnCampbell*cos(sza)+EdhCampbell
	
    ###Quality control
    lower=0
    EbnCampbell[EbnCampbell<lower]=0
    EdhCampbell[EdhCampbell<lower]=0
    EghCampbell[EghCampbell<lower]=0
	return(list(EbnCampbell, EdhCampbell, EghCampbell))
	}


