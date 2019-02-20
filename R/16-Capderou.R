#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 16-Capderou 1987
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Capderou, M. (1987). Theoretical and experimental models solar atlas of Algeria (in French) Tome 1 and 2. Algeria: University Publications Office.
#Marif, Y., Chiba, Y., Belhadj, M. M., Zerrouki, M., & Benhammou, M. (2018). A clear sky irradiation assessment using a modified Algerian solar atlas model in Adrar city. Energy Reports, 4, 84-90.

###Inputs:
#  Esc=1367   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  latitude   [double]    (location's latitude，-90~90)
#  altitude   [m]         (location's altitude)
#  press      [mb]        (local barometric)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceCapderou<-function(latitude,altitude){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
	B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #air mass
    mA=press/1013.25*(cos(sza)+0.15*((3.885+(pi/2-sza)*180/pi)^(-1.253)))^(-1)
    mAc=(cos(sza)+0.15*((3.885+(pi/2-sza)*180/pi)^(-1.253)))^(-1)
    
	Ahe=sin(360/365*(Dayth-121))
	#the atmospheric turbidity caused by water vapor absorption 
	T0=2.4-0.9*sin(latitude)+0.1*Ahe*(2+sin(latitude))-0.2*altitude/1000-(1.22+0.14*Ahe)*(1-cos(sza))
	#the atmospheric turbidity corresponding to the molecular diffusion 
	T1=0.89^(altitude/1000)
	#the atmospheric turbidity relative to the aerosol diffusion coupled with a slight absorption 
	T2=(0.9+0.4*Ahe)*0.63^(altitude/1000)
	#the atmospheric Linke turbidity factor 
	Tlc=T0+T1+T2
	
	delta_rk=1/(9.4+0.9*mA)
	#direct normal irradiance
	EbnCapderou=Eext*exp(-Tlc*mAc*delta_rk)
	
	#diffuse horizontal irradiance
	a=1.1
	b=log(T1+T2)-2.8+1.02*(1-cos(sza))^2
	EdhCapderou=Eext*exp(-1+1.06*log(cos(sza))+a-sqrt(a^2+b^2))
	#global horizontal irradiance
	EghCapderou=EbnCapderou*cos(sza)+EdhCapderou
	
    ###Quality control
    lower=0
    EbnCapderou[EbnCapderou<lower]=0
    EdhCapderou[EdhCapderou<lower]=0
    EghCapderou[EghCapderou<lower]=0
	return(list(EbnCapderou, EdhCapderou, EghCapderou))
	}

