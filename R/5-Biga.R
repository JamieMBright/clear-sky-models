#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 05-Biga 1979
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Schulze, R. E. (1976). A physically based method of estimating solar radiation from suncards. Agricultural Meteorology, 16(1), 85-101.

###Inputs:
#  sza  [radians]  (zenith_angle) 

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Codes:
IrradianceBiga<-function(){
	#direct normal irradiance
	EbnBiga=926*((cos(sza))^0.29)
	#diffuse horizontal irradiance
	EdhBiga=131*((cos(sza))^0.6)
	#global horizontal irradiance
	EghBiga=1057*((cos(sza))^1.17)
	
    ###Quality control
    lower=0
    EbnBiga[EbnBiga<lower]=0
    EdhBiga[EdhBiga<lower]=0
    EghBiga[EghBiga<lower]=0
	return(list(EbnBiga, EdhBiga, EghBiga))
	}


