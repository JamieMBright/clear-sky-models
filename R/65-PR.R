#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 65-PR 2000
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Gueymard, C. A. (2003). Direct solar transmittance and irradiance predictions with broadband models. Part I: detailed theoretical performance assessment. Solar Energy, 74(5), 355-379.
#Psiloglou, B. E., Santamouris, M., & Asimakopoulos, D. N. (2000). Atmospheric broadband model for computation of solar radiation at the Earth’s surface. Application to Mediterranean climate. pure and applied geophysics, 157(5), 829-860.
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

###Inputs:
#  Esc=1373     [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  press        [mb]               (local barometric)
#  albedo       [double]           (surface/ground albedo)
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
IrradiancePR<-function(){	
	#extraterrestrial irradiance
	Esc=1373
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))	
	
    #air mass
    amky=(cos(sza)+0.50572*((pi/2-sza)/pi*180+6.07995)^-1.6364)^-1
    amkyp=amky*press/1013.25
    
    #transmittances from Psiloglou2000
    #Badescu added the presure correction where needed for cinsistency
    Trr=exp(-0.1128*(amkyp^0.8346)*(0.9341-amkyp^0.9868+0.9391*amkyp)) 
    ako3=amky*ozone
    Tro=1-0.2554*ako3/((1+6107.26*ako3)^0.204+0.471*ako3)
    akco2=amkyp*330
    TrCO2=1-0.0721*akco2/((1+377.89*akco2)^0.5855+3.1709*akco2)
    akco=amkyp*0.075
    TrCO=1-0.0062*akco/((1+243.67*akco)^0.4246+1.7222*akco)
    akch4=amkyp*1.6
    TrCH4=1-0.0192*akch4/((1+166.095*akch4)^0.4221+0.7186*akch4)
    akn2o=amkyp*0.28
    TrN2O=1-0.0326*akn2o/((1+107.413*akn2o)^0.5501+0.9093*akn2o)
    ako2=amkyp*2.095e5
    TrO2=1-0.0003*ako2/((1+476.934*ako2)^0.4892+0.1261*ako2)
    Trmg=TrCO2*TrCO*TrCH4*TrN2O*TrO2
    akw=amky*wv
    Trw=1-3.014*akw/((1+119.3*akw)^0.644+5.814*akw) 
    
    #Aerosol transmittance from REST 
    Masa=amky
    beta2=ang_beta*ang_beta
    b0=1.6933
    bx=1+0.42003*ang_beta
    b1=(-0.013029+0.13126*ang_beta)/bx
    b2=(-0.0083581+0.40323*ang_beta+0.123*beta2)/bx
    TauA=ang_beta*(b0+Masa*b1)/(1+b2*Masa)
    Tra=exp(-Masa*TauA) 

	#direct normal radiation
	EbnPR=Eext*Trr*Trmg*Tro*Trw*Tra
	Eb=EbnPR*cos(sza)	
	
	#diffuse horizontal radiation
	Traa=1-1.405e-3*amky-9.013e-5*amky^2+2.2e-6*amky^3
	Tras=Tra/Traa
	ros=0.08503
	Ed1=Eext*cos(sza)*Trmg*Tro*Trw*Traa*(1-Tras*Trr)/2
	Ed2=(Ed1+Eb)*albedo*ros/(1-albedo*ros)
	EdhPR=Ed1+Ed2 
	
	#global horizontal irradiance
	EghPR=EbnPR*cos(sza)+EdhPR
	
	###Quality control
    lower=0
    EbnPR[EbnPR<lower]=0
    EdhPR[EdhPR<lower]=0
    EghPR[EghPR<lower]=0
	return(list(EbnPR, EdhPR, EghPR))
	}


