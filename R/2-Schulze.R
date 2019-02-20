#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 02-Schulze 1957
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
IrradianceSchulze<-function(){
	#direct normal irradiance
	Ebnschulze=1127*(0.888^(1/cos(sza)))
	
	#diffuse horizontal irradiance
	Edhschulze=94.23*((cos(sza))^0.5)
	
	#global horizontal irradiance
	Eghschulze=Ebnschulze*cos(sza)+Edhschulze
	
    ###Quality control
    lower=0
    Ebnschulze[Ebnschulze<lower]=0
    Edhschulze[Edhschulze<lower]=0
    Eghschulze[Eghschulze<lower]=0
	return(list(Ebnschulze, Edhschulze, Eghschulze))
	}
	
	