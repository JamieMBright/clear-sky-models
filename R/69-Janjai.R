#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 69-Janjai 2011
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Janjai, S., Sricharoen, K., & Pattarapanitchai, S. (2011). Semi-empirical models for the estimation of clear sky solar global and direct normal irradiances in the tropics. Applied Energy, 88(12), 4749-4755.

###Inputs:
#  Esc=1366.1   [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  altitude     [m]                (location's altitude)
#  ang_alpha    [dimensionless]    (Angstrom_exponent, also known as alpha)
#  ang_beta     [dimensionless]    (Angstrom turbidity coefficient)
#  ozone        [atm.cm]           (total columnar amount of ozone)
#  wv           [atm.cm]           (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceJanjai<-function(altitude){
    #extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #air mass
    am=1/(cos(sza)+0.15*((pi/2-sza)/pi*180+3.885)^(-1.253))
    ama=am*exp(-0.0001184*altitude)
    
    #efficients
    a1 = 0.778227
    a2 = 0.71640 
    b1 = 1.198932
    b2 = 0.35320 
    c1 = -0.106634
    c2 = 0.10126 
    d1 = 0.337373
    d2 = 0.841372 
    e1 = 0.009181
    e2 = 0.017649 
    f1 = -0.009852
    f2 = 0.004851 
    g1 = 0.482012
    g2 = -0.48286
    
   #direct norm irradiance
   A2=a2*Eext*(cos(sza))^b2
   B2=c2+d2*ang_beta+e2*ang_alpha+f2*wv+g2*ozone
   Ebnjanjai=A2*exp(-B2*ama)
   
    #global horizontal irradiance 
    A1=a1*Eext*(cos(sza))^b1
    B1=c1+d1*ang_beta+e1*ang_alpha+f1*wv+g1*ozone
    Eghjanjai=A1*exp(-B1*ama)
    
    #diffuse horizontal irradiance
    Edhjanjai=Eghjanjai-Ebnjanjai*cos(sza)
    
	###Quality control
    lower=0
    Ebnjanjai[Ebnjanjai<lower]=0
    Edhjanjai[Edhjanjai<lower]=0
    Eghjanjai[Eghjanjai<lower]=0
	return(list(Ebnjanjai,Edhjanjai,Eghjanjai))
	}
	
	
	