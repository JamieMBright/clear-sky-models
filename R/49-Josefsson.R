#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 49-Josefsson 1985
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  Esc=1366.1   [Wm-2]      (Solar constant)
#  sza          [radians]   (zenith_angle) 
#  press        [mb]        (local barometric)
#  wv           [atm.cm]    (total columnar amount of water vapour)
#  albedo       [double]    (surface/ground albedo)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceJosefsson<-function(){	
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B)	)
	
	#air mass
	amk=1/(cos(sza)+0.15*((pi/2-sza)/pi*180+3.885)^(-1.253))
    am=amk*press/1013.25
    am2=am*am
    am3=am2*am
    am4=am3*am 
    am5=am4*am 
    am66=1.66*press/1013.25 
    
    h=(pi/2-sza)*180/pi
    w1=10*wv 
    omeg=0.75
    wp=w1*amk
    aw=0.29*wp/((1+14.15*wp)^0.635+0.5925*wp)
    g=0.5248+0.007981*h
    h_xiabiao=which(h>45)
    g[h_xiabiao]=0.856+0.000734*h[h_xiabiao]
    
    #transmittances
    Tra=0.95^amk 
    Trab=0.95^am66
    Traa=1-(1-omeg)*(1-Tra)
    Tras=1-omeg*(1-Tra)
    TrO3=0.95545
    TrR=0.9768-0.0874*am+0.010607552*am2-0.000846205*am3+3.57246e-5*am4-6.0176e-7*am5
    Tabs=TrO3*Traa	

	#direct normal radiation
	EbnJosefsson=Eext*(Tabs*Tras*TrR-aw)
	Eb=EbnJosefsson*cos(sza)
	
	#diffuse horizontal radiation
	Fa=0.93-0.21*log(amk)
	Fab=0.93-0.21*log(1.66)
	rhos=0.0685+(1-Trab)*omeg*(1-Fab)
	Ed0=Eext*cos(sza)*(0.5*Tabs*Tras*(1-TrR)+g*(TrR*Tabs-aw)*(1-Tras)) 
	Edb=albedo*rhos*(Eb+Ed0)/(1-albedo*rhos) 
	EdhJosefsson=Ed0+Edb	
	
	#global radiation
	EghJosefsson=EbnJosefsson*cos(sza)+EdhJosefsson
	
	###Quality control
    lower=0
    EbnJosefsson[EbnJosefsson<lower]=0
    EdhJosefsson[EdhJosefsson<lower]=0
    EghJosefsson[EghJosefsson<lower]=0
	return(list(EbnJosefsson, EdhJosefsson, EghJosefsson))
	}
	
	
	