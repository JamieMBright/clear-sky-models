#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 09-Yang & Walsh  2014
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Yang, D., Walsh, W. M., & Jirutitijaroen, P. (2014). Estimation and applications of clear sky global horizontal irradiance at the equator. Journal of Solar Energy Engineering, 136(3), 034505.

###Inputs:
#  Esc=1362   [Wm-2]     (Solar constant)
#  sza        [radians]  (zenith_angle) 

###Outputs:
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceYangWalsh<-function(){
	#extraterrestrial irradiance
	Esc=1362
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
	#global horizontal irradiance
	Eghyangwalsh=0.8298*Eext*((cos(sza))^1.3585)*exp(-0.00135*(pi/2-sza)/pi*180)
		
    ###Quality control
    lower=0
    Eghschulze[Eghyangwalsh<lower]=0
	return(Eghyangwalsh)
	}
	
