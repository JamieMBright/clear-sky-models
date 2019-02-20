#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 47-CLS 1976
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Suckling, P. W., & Hay, J. E. (1976). Modelling direct, diffuse, and total solar radiation for cloudless days. Atmosphere, 14(4), 298-308.
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  Esc=1353   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  press      [mb]        (local barometric)
#  wv         [atm.cm]    (total columnar amount of water vapour)
#  albedo     [double]    (surface/ground albedo)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceCLS<-function(){	
	#extraterrestrial irradiance
	Esc=1353
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))	
	
	#air mass
	amk=1/(cos(sza)+0.15*((pi/2-sza)/pi*180+3.885)^(-1.253))
    am=amk*press/1013.25
	am2=am*am
	am3=am2*am
	am4=am3*am
	
	#water vapour
	wm=wv*am
	wmb=wv*1.66 
	
	#transmittances
	TrR=0.972-0.08262*am+0.00933*am2-0.00095*am3+4.37e-5*am4
	Trrb=0.85655
	Tws=1-0.0225*wm
	Twsb=1-0.0225*wmb
	Twa=1-0.077*(wm^0.3)
	Twab=1-0.077*(wmb^0.3)
	Tda=0.975^am
	Tdab=0.975^1.66
	Tds=Tda
	Tdsb=Tdab
	
	#direct normal radiation	
	EbnCLS=Eext*TrR*Tda*Tds*Twa*Tws
	
	#diffuse horizontal radiation
	Eds=0.6*Eext*cos(sza)*Twa*Tda*(1-Tws*Tds*TrR)
	Es=EbnCLS*cos(sza)+Eds
	Edb=albedo*Es*0.4*Twab*Tdab*(1-Twsb*Tdsb*Trrb)
	EdhCLS=Eds+Edb
	
	#global radiation
	EghCLS=Es+Edb
	
	###Quality control
    lower=0
    EbnCLS[EbnCLS<lower]=0
    EdhCLS[EdhCLS<lower]=0
    EghCLS[EghCLS<lower]=0
	return(list(EbnCLS, EdhCLS, EghCLS))
	}
	
	