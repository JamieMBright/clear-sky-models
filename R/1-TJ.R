#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 01-TJ 1957
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Threlkeld, J. L., & Jordan, R. C. (1957). Direct solar radiation available on clear days. Heat., Piping Air Cond., 29(12).
#Masters, G. M. (2013). Renewable and efficient electric power systems. John Wiley & Sons.

###Inputs:
#  sza  [radians]  (zenith_angle) 

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceTJ<-function(){  
    #air mass
    mA=1/cos(sza)
    
    #coefficients
    A=1160+75*sin(360*(Dayth-275)/365)
    k=0.174+0.035*sin(360*(Dayth-100)/365)
	C=0.095+0.04*sin(360*(Dayth-100)/365)
	
	#direct normal irradiance
	EbnTJ=A*exp(-k*mA)
	
	#diffuse horizontal irradiance
	EdhTJ=C*EbnTJ
	
	#global horizontal irradiance
	EghTJ=EbnTJ*cos(sza)+EdhTJ
	
    #Quality control
    lower=0
    EbnTJ[EbnTJ<lower]=0
    EdhTJ[EdhTJ<lower]=0
    EghTJ[EghTJ<lower]=0
	return(list(EbnTJ, EdhTJ, EghTJ))
	}

