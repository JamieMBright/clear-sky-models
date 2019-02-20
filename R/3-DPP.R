#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 03-DPP 1978
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Daneshyar, M. (1978). Solar radiation statistics for Iran. Sol. Energy;(United States), 21(4).
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  sza  [radians]  (zenith_angle) 

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Codes:
IrradianceDPP<-function(){
	#direct normal irradiance
	Ebndpp=950*(1-exp(-0.075*(pi/2-sza)/pi*180))
	
	#diffuse horizontal irradiance 
	Edhdpp=2.534+3.475*(pi/2-sza)/pi*180
	
	#global horizontal irradiance
	Eghdpp=Ebndpp*cos(sza)+Edhdpp
	
	###Quality control
    lower=0
    Ebndpp[Ebndpp<lower]=0
    Edhdpp[Edhdpp<lower]=0
    Eghdpp[Eghdpp<lower]=0
	return(list(Ebndpp, Edhdpp, Eghdpp))	
	}

