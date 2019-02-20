#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 15-KASM  1983 
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Kasten, F. (1984). Parametrisierung der Globalstahlung durch Bedeckungsgrad und Trubungsfaktor. Annalen der Meteorologie Neue, 20, 49-50.
#Badescu, V. (1997). Verification of some very simple clear and cloudy sky models to evaluate global solar irradiance. Solar Energy, 61(4), 251-264.
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  Esc=1367   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  press      [mb]        (local barometric)
#  wv         [atm.cm]    (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceKASM<-function(){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
	B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #air mass
    amk=(cos(sza)+0.15*((3.885+(pi/2-sza)*180/pi)^(-1.253)))^(-1)
    am=amk*press/1013.25
    am2=am*am
    am3=am2*am
    am4=am3*am
    am5=am4*am
    
    w1=10*wv
    TrR=0.9768-0.0874*am+0.010607552*am2-0.000846205*am3+3.57246e-5*am4-6.0176e-7*am5
    
    x=amk*3.5 
    aov=0.002118*x/(1+0.0042*x+0.00000323*x*x)
    aou=0.1082*x/((1+13.86*x)^0.805)+0.00658*x/(1+(10.36*x)^3)
    TroP=1-aov-aou
    
    wp=w1*am
    aw=0.29*wp/((1+14.15*wp)^0.635+0.5925*wp)
    Tra=0.9^am 
    
    #direct normal irradiance
	EbnKASM=Eext*(TrR*TroP-aw)*Tra
	
	#diffuse horizontal irradiance
	EdhKASM=Eext*cos(sza)*(TroP*(1-TrR)*0.5+(TroP*TrR-aw)*(1-Tra)*0.75*(0.93-0.21*log(am))) 

	#global horizontal irradiance
	EghKASM=EbnKASM*cos(sza)+EdhKASM
	
    ###Quality control
    lower=0
    EbnKASM[EbnKASM<lower]=0
    EdhKASM[EdhKASM<lower]=0
    EghKASM[EghKASM<lower]=0
	return(list(EbnKASM, EdhKASM, EghKASM))
	}
	
	