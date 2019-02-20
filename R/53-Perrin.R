#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 53-Perrin 1975
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#de Brichambaut, C. P. (1975). Estimation des Ressources Energétiques en France. Cahiers de l’AFEDES, (1).
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  Esc=1366.1   [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  press        [mb]               (local barometric)
#  ang_beta     [dimensionless]    (Angstrom turbidity coefficient)
#  wv           [atm.cm]           (total columnar amount of water vapour)
#  ozone        [atm.cm]           (total columnar amount of ozone)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradiancePerrin<-function(){
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))  
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #air mass
	am=1/(cos(sza)+0.15*((pi/2-sza)/pi*180+3.885)^(-1.253))
    ama=am*press/1013.25
    
    Trr=exp(-0.031411-0.064331*ama)
    
    ako3=am*ozone 
    abso3=0.015+0.024*ako3    
    Tro=1-abso3
    	
	akw=am*wv
    xw=vector(mode='numeric',length=length(wv)) #create a zero vector
    xiabiao_of_xw=which(akw>0)
    xw[xiabiao_of_xw]=log(akw[xiabiao_of_xw])
    #xw[akw>0]=log(akw[akw>0])
    xw2=xw*xw
    absw=0.1+0.03*xw+0.002*xw2
    Trw=1-absw
    
    xwg=vector(mode='numeric',length=length(wv)) #create a zero vector    
    akwg=ama*wv
    xiabiao_of_xwg=which(akwg>0)
    xwg[xiabiao_of_xwg]=log(akwg[xiabiao_of_xwg])
    #xwg[akwg>0]=log(akwg[akwg>0])
    absg=0.013-0.0015*xwg
    Trg=1-absg
    
    Tra=exp(-1.4327*am*ang_beta)

	#direct normal irradiance
	EbnPerrin=Eext*Trr*Tra*(1-abso3-absw-absg)
	
	tCDA=1/(9.4+0.9*ama) 
	TL=-(log(EbnPerrin/Eext))/(tCDA*ama)
	#global horizontal irradiance
	EghPerrin=(1270-56*TL)*(cos(sza))^((TL+36)/33)
	
	#diffuse horizontal irradiance
	EdhPerrin=EghPerrin-EbnPerrin*cos(sza)
	
	#Strange limit in Badescu
	xiabiao_of_egh=which(EdhPerrin<0)
    EghPerrin[xiabiao_of_egh]=EbnPerrin[xiabiao_of_egh]*(cos(sza))[xiabiao_of_egh]
    #EghPerrin[EdhPerrin<0]=EbnPerrin[EdhPerrin<0]*(cos(sza))[EdhPerrin<0]

    ###Quality control
    lower=0
    EbnPerrin[EbnPerrin<lower]=0
    EdhPerrin[EdhPerrin<lower]=0
    EghPerrin[EghPerrin<lower]=0
	return(list(EbnPerrin, EdhPerrin, EghPerrin))
	}

