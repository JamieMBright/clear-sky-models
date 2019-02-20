#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 04-Adnot 1979
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Adnot, J., Bourges, B., Campana, D., & Gicquel, R. (1979). Utilisation de courbes de fréquences cumulées d'irradiation solaire globale pour le calcul des installations solaires.
#Barbieri, F., Rifflart, C., Vo, B. T., Rajakaruna, S., & Ghosh, A. (2016, July). A comparative study of clear-sky irradiance models for Western Australia. In 2016 IEEE Power and Energy Society General Meeting (PESGM) (pp. 1-5). IEEE.

###Inputs:
#  sza  [radians]  (zenith_angle) 

###Outputs:
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Codes:
IrradianceAdnot<-function(){
	#global horizontal irradiance
	EghAdnot=951.39*(cos(sza))^1.15
	
    ###Quality control
    lower=0
    EghAdnot[EghAdnot<lower]=0
	return(EghAdnot)
	}

