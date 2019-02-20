#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 48-King 1979
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#King, R., & Buckius, R. O. (1979). Direct solar transmittance for a clear sky. Solar Energy, 22(3), 297-301.
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  Esc=1366.1   [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  press        [mb]               (local barometric)
#  albedo       [double]           (surface/ground albedo)
#  wv           [atm.cm]           (total columnar amount of water vapour)
#  ang_beta     [dimensionless]    (Angstrom turbidity coefficient)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceKing<-function(){	
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))	
    	
	#air mass
	amk=1/(cos(sza)+0.15*((pi/2-sza)/pi*180+3.885)^(-1.253))
    am=amk*press/1013.25
    amb=amk*ang_beta 
    a1=0.5
    
    Tscat=(5.228/am)*(exp(-0.0002254*am)-exp(-0.1409*am))*exp(-amb/(0.6489^1.3))+0.00165*am+0.2022*(1-amb/(2.875^1.3))
    
    Tabs=(0.1055+0.07053*log10(amk*wv+0.07854))*exp(-amb/(1.519^1.3))
 	
 	#direct normal radiation   
    EbnKing=Eext*(Tscat-Tabs)
    
    taul=0.8336+0.17257*ang_beta-0.64595*exp(-7.4296*ang_beta^1.5)
    Tml=exp(-am*taul)
    F1=2*(1+2*cos(sza)+(1-2*cos(sza))*Tml)
    F2=(albedo-1)*(a1*taul-4*taul-4)+4*albedo
    Tdif=0.634*(F1/F2-Tml)
    #diffuse horizontal radiation
    EdhKing=Eext*cos(sza)*Tdif
	
	#global radiation
	EghKing=EbnKing*cos(sza)+EdhKing
	
	###Quality control
    lower=0
    EbnKing[EbnKing<lower]=0
    EdhKing[EdhKing<lower]=0
    EghKing[EghKing<lower]=0
	return(list(EbnKing, EdhKing, EghKing))
	}


