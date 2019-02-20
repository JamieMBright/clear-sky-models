#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 67-MRMv5 2008
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Kambezidis, H. D., & Psiloglou, B. E. (2008). The meteorological radiation model (MRM): advancements and applications. In Modeling solar radiation at the earth’s surface (pp. 357-392). Springer, Berlin, Heidelberg.

###Inputs:
#  Esc=1366.1     [Wm-2]             (Solar constant)
#  sza            [radians]          (zenith_angle) 
#  press          [mb]               (local barometric)
#  albedo         [double]           (surface/ground albedo)
#  ang_beta       [dimensionless]    (Angstrom turbidity coefficient)
#  ozone          [atm.cm]           (total columnar amount of ozone)
#  wv             [atm.cm]           (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceMRM5<-function(){
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
        
    #air mass
    am=(cos(sza)+0.50572*((pi/2-sza)/pi*180+6.07995)^-1.6364)^-1
    ama=am*press/1013.25

    #atmospheric gases transmittance function
    Tco2=1-( 0.07210*ama*350/((1+377.890*ama*350)^0.5855+3.1709*ama*350) )
    Tco=1-( 0.0062*ama*0.075/((1+243.670*ama*0.075)^0.4246+1.7222*ama*0.075) )
    Tn2o=1-( 0.0326*ama*0.28/((1+107.413*ama*0.28)^0.5501+0.9093*ama*0.28) )
    Tch4=1-( 0.0192*ama*1.6/((1+166.095*ama*1.6)^0.4221+0.7186*ama*1.6) )
    To2=1-( 0.0003*ama*2.095e5/((1+476.934*ama*2.095e5)^0.4892+0.1261*ama*2.095e5) )
    #the transmittance function for the mixed gases
    Tmg=Tco2*Tco*Tn2o*Tch4*To2
    
    #the transmittance function for water vapour  #because we have wv, use eq.14.27,otherwise eq.14.20abcd
    lw=wv
    #Tdok=(Tem+273.15)/100
    #ess=exp(22.329699-49.140396*Tdok^-1-10.921853*Tdok^-2-0.39015156*Tdok)
    #lw=0.00493*ess*RH/Tdok
    Tw=1-( 3.0140*am*lw/((1+119.300*am*lw)^0.6440+5.8140*am*lw) )
    
    #the transmittance function for ozone  #because we have ozone, use eq.14.27,otherwise eq.14.9
    To=1-( 0.2554*am*ozone/((1+6107.260*am*ozone)^0.2040+0.4710*am*ozone) )
    
    #the transmittance function for Rayleigh scattering    
    Tr=exp(-0.1128*ama^0.8346*(0.9341-ama^0.9868+0.9391*ama))
    
    #the Mie scattering transmittance function
    Ta=exp(-am*ang_beta*(0.6777+0.1464*(am*ang_beta)-0.00626*(am*ang_beta)^2)^-1.3)
    
    #the aerosol absorption function
    Taa=1-0.1*(1-am+am^1.06)*(1-Ta)
    
    #the aerosol scattering function
    Tas=Ta/Taa
    
    #direct beam irradiance
    EbnMRM5=Eext*Ta*Tr*To*Tw*Tmg
    
    #the circumsolar diffuse radiation produced by a single-scattering mode of molecules and aerosols 
    Ids=Eext*cos(sza)*Taa*To*Tw*Tmg*0.5*(1-Tas*Tr)     
    #the diffuse component reflected by the ground and backscattered by the atmosphere        
    Ta166=exp(-1.66*ang_beta*(0.6777+0.1464*(1.66*ang_beta)-0.00626*(1.66*ang_beta)^2)^-1.3)
    pa=0.0685+0.16*(1-Ta166)   #sky reflectance
    Idm=(EbnMRM5*cos(sza)+Ids)*(albedo*pa/(1-albedo*pa))
    #diffuse horizontal irradiance
    EdhMRM5=Ids+Idm
    
    #global horizontal irradiance
    EghMRM5=EbnMRM5*cos(sza)+EdhMRM5
    
	###Quality control
    lower=0
    EbnMRM5[EbnMRM5<lower]=0
    EdhMRM5[EdhMRM5<lower]=0
    EghMRM5[EghMRM5<lower]=0
	return(list(EbnMRM5, EdhMRM5, EghMRM5))
	}
	
	
	
