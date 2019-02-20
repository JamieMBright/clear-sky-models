#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 59-MMAC 1993,2003
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Gueymard, C. (1993). Critical analysis and performance assessment of clear sky solar irradiance models using theoretical and measured data. Solar Energy, 51(2), 121-138.
#Gueymard, C. A. (2003). Direct solar transmittance and irradiance predictions with broadband models. Part I: detailed theoretical performance assessment. Solar Energy, 74(5), 355-379.
#Davies, J. A., & McKay, D. C. (1982). Estimating solar irradiance and components. Solar Energy, 29(1), 55-64.

###Inputs:
#  Esc=1353         [Wm-2]              (Solar constant)
#  sza              [radians]           (zenith_angle) 
#  press            [mb]                (local barometric)
#  albedo           [double]            (surface/ground albedo)
#  Broadband_AOD    [dimensionless]     (broadband aerosol optical depth)
#  wv               [atm.cm]            (total columnar amount of water vapour)


###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceMMac<-function(){
	#extraterrestrial irradiance
	Esc=1353
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))

	#Air Mass
    m=35/((1224*(cos(sza))^2+1)^0.5)
    m[m<0.0]=NA 
    ma=m*press/1013.25
    
    #aerosol transmittance
    Ta=exp(-ma*broadbandaod)
    
    #Rayleigh transmittance
    TR=1/(1+1/exp(2.12182-0.791532*log(ma)+0.024761*(log(ma))^2))

    #Ozone Transmittance   
    XO=m*(0.35*10) #Davies and Mckay 1982 set ozone a fixed value of 3.5mm,here in the code unit is cm    
    AO=((0.1082*XO)/((1+13.86*XO)^0.805))+((0.00658*XO)/(1+(10.36*XO)^3))+(0.002118*XO/(1+0.0042*XO+3.23e-6*XO^2))
    TO=1-AO
    
    #Water vapor transmittance
    XW=10*m*wv*(press/1013.25)^0.75
    AW=0.29*XW/((1+ 14.15*XW)^0.635+0.5925*XW)
    
    #direct beam irradiance
    Ebnmmac=Eext*(TO*TR-AW)*Ta
    
    #Forward aerosol scatterance
    Ba=0.93-0.21*log(m)
    
    #Aerosol single-scattering albedo
    wa=0.98
    Taaa=0.95^1.66
    Baa=0.93-0.21*log(1.66)
    #diffuse components from scattering by aerosol
    Ida=Eext*cos(sza)*wa*Ba*(TO*TR-AW)*(1-Ta)    
    #diffuse components from Rayleigh scatter
    Idr=Eext*cos(sza)*(0.5*TO*Taaa*(1-TR))
    #diffuse horizontal irradiance without albedos
    Edh=Idr+Ida
    #diffuse horizontal irradiance with albedos
    poub=0.0685+(1-Taaa)*wa*(1-Baa)
    Edhmmac=poub*albedo*(Ebnmmac*cos(sza)+Edh)/(1-poub*albedo)+Edh
    
    #global horizontal irradiance
    Eghmmac=(Ebnmmac*cos(sza)+Edh)/(1-poub*albedo)

    ###Quality control
    lower=0
    Ebnmmac[Ebnmmac<lower]=0
    Edhmmac[Edhmmac<lower]=0
    Eghmmac[Eghmmac<lower]=0
	return(list(Ebnmmac, Edhmmac, Eghmmac))
	}

