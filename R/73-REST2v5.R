#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 73-REST2  version5  2008
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Gueymard, C. A. (2008). REST2: High-performance solar radiation model for cloudless-sky irradiance, illuminance, and photosynthetically active radiation–Validation with a benchmark dataset. Solar Energy, 82(3), 272-285.

###Inputs:
#  Esc=1366.1   [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  press        [mb]               (local barometric)
#  albedo       [double]           (surface/ground albedo)
#  ang_alpha    [dimensionless]    (Angstrom_exponent, also known as alpha)
#  ang_beta     [dimensionless]    (Angstrom turbidity coefficient)
#  ozone        [atm.cm]           (total columnar amount of ozone)
#  NO2          [atm.cm]           (total columnar amount of nitrogen dioxide)
#  wv           [atm.cm]           (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceRest2v5<-function(){
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))  
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #air mass for aerosols extinction 
    #Gueymard2013<Direct solar transmittance and irradiance predictions with broadband models. Part I: detailed theoretical performance assessment > appendix B.8-B.12
    ama=(cos(sza)+0.16851*(sza*180/pi)^0.18198/(95.318-sza*180/pi)^1.9542)^-1
    #air mass for water vapor absorption
    amw=(cos(sza)+0.10648*(sza*180/pi)^0.11423/(93.781-sza*180/pi)^1.9203)^-1
    #air mass for nitrogen dioxide absorption
    amn=(cos(sza)+1.1212*(sza*180/pi)^1.6132/(111.55-sza*180/pi)^3.2629)^-1
    #air mass for ozone absorption
    amo=(cos(sza)+1.0651*(sza*180/pi)^0.6379/(101.8-sza*180/pi)^2.2694)^-1
    #air mass for Rayleigh scattering and uniformly mixed gases absorption
    amR=(cos(sza)+0.48353*((sza*180/pi)^0.095846)/(96.741-sza*180/pi)^1.754)^-1
    amRe=(press/1013.25)*(cos(sza)+0.48353*((sza*180/pi)^0.095846)/(96.741-sza*180/pi)^1.754)^-1
  
    ## Band 1 ##
    #transmittance for Rayleigh scattering 
    TR1=(1 + 1.8169*amRe - 0.033454*amRe^2 )/(1 + 2.063*amRe + 0.31978*amRe^2)  
    #transmittance for uniformly mixed gases absorption
    Tg1=(1 + 0.95885*amRe + 0.012871*amRe^2) /(1 + 0.96321*amRe + 0.015455*amRe^2)
    #transmittance for Ozone absorption
    uo=ozone
    f1=uo*(10.979 - 8.5421*uo)/(1 + 2.0115*uo + 40.189*uo^2)
    f2=uo*(-0.027589 - 0.005138*uo)/(1 - 2.4857*uo + 13.942*uo^2)
    f3=uo*(10.995 - 5.5001*uo)/(1 + 1.6784*uo + 42.406*uo^2)
    To1=(1 + f1*amo + f2*amo^2)/(1 + f3*amo)     
    #transmittance for Nitrogen dioxide absorption
    un=NO2
    g1=(0.17499 + 41.654*un - 2146.4*un^2)/(1 + 22295*un^2)
    g2=un*(-1.2134 + 59.324*un)/(1 + 8847.8*un^2)
    g3=(0.17499 + 61.658*un + 9196.4*un^2)/(1 + 74109*un^2)
    Tn1=pmin(1,((1 + g1*amw +  g2*amw^2)/(1 + g3*amw)))
    Tn1166=pmin(1,((1 + g1*1.66 +  g2*1.66^2)/(1 + g3*1.66)))   #at air mass=1.66
    #transmittance for Water Vapor absorption
    h1=wv*(0.065445 + 0.00029901*wv)/(1 + 1.2728*wv) 
    h2=wv*(0.065687 + 0.0013218*wv)/(1 + 1.2008*wv)
    Tw1=(1 + h1*amw)/(1 + h2*amw)
    Tw1166=(1 + h1*1.66)/(1 + h2*1.66)     #at air mass=1.66

     #coefficients of angstrom_alpha
     AB1=ang_beta
     alph1=ang_alpha
     d0=0.57664 - 0.024743*alph1
     d1=(0.093942 - 0.2269*alph1 + 0.12848*alph1^2)/(1 + 0.6418*alph1)
     d2=(-0.093819 + 0.36668*alph1 - 0.12775*alph1^2)/(1 - 0.11651*alph1)
     d3=alph1*(0.15232 - 0.087214*alph1 + 0.012664*alph1^2)/(1 - 0.90454*alph1 + 0.26167*alph1^2)
     ua1=log(1 + ama*AB1)   #below Eq.7a
     lam1=(d0 + d1*ua1 + d2*ua1^2)/(1 + d3*ua1^2)
     
     #Aeroso transmittance   Eq.6,7
     ta1=AB1*lam1^-alph1
     TA1=exp(-ama*ta1)
     #Aeroso scattering transmittance
     TAS1=exp(-ama*0.92*ta1)    #w1=0.92 recommended by author

     #forward scattering fractions for Rayleigh extinction    Eq.10,11
     BR1=0.5*(0.89013 - 0.0049558*amR + 0.000045721*amR^2)
     #the aerosol forward scatterance factor
     Ba=1 - exp(-0.6931-1.8326*cos(sza))

     #Aerosol scattering correction factor      Appendix
     g0=(3.715 + 0.368*ama + 0.036294*ama^2)/(1 + 0.0009391*ama^2)
     g1=(-0.164 - 0.72567*ama + 0.20701*ama^2)/(1 + 0.0019012*ama^2)
     g2=(-0.052288 + 0.31902*ama + 0.17871*ama^2)/(1 + 0.0069592*ama^2)
     F1=(g0 + g1*ta1)/(1 + g2*ta1)
     
     #sky albedo     Appendix
     rs1=(0.13363 + 0.00077358*alph1 + AB1*(0.37567 + 0.22946*alph1)/(1-0.10832*alph1))/(1 + AB1*(0.84057 + 0.68683*alph1)/(1 - 0.08158*alph1)) 
     #ground albedo
     rg=albedo

     ## Band 2 ##
     #transmittance for Rayleigh scattering 
     TR2=(1 - 0.010394*amRe)/(1-0.00011042*amRe^2)
     #transmittance for uniformly mixed gases absorption
     Tg2=(1 + 0.27284*amRe - 0.00063699*amRe^2)/(1 + 0.30306*amRe)
     #transmittance for Ozone absorption
     To2=1 # Ozone (none)
     #transmittance for Nitrogen dioxide absorption
     Tn2=1 # Nitrogen (none)
     Tn2166=1     #at air mass=1.66
     #transmittance for water vapor  absorption
     c1=wv*(19.566 - 1.6506*wv + 1.0672*wv^2)/(1 + 5.4248*wv + 1.6005*wv^2)
     c2=wv*(0.50158 - 0.14732*wv + 0.047584*wv^2)/(1 + 1.1811*wv + 1.0699*wv^2)
     c3=wv*(21.286 - 0.39232*wv + 1.2692*wv^2)/(1 + 4.8318*wv + 1.412*wv^2)
     c4=wv*(0.70992 - 0.23155*wv + 0.096514*wv^2)/(1 + 0.44907*wv + 0.75425*wv^2)
     Tw2=(1 + c1*amw + c2*amw^2)/(1 + c3*amw + c4*amw^2)
     Tw2166=(1 + c1*1.66 + c2*1.66^2)/(1 + c3*1.66 + c4*1.66^2)

     #coefficients of angstrom_alpha
     AB2=ang_beta
     alph2=ang_alpha   
     e0=(1.183 - 0.022989*alph2 + 0.020829*alph2^2)/(1 + 0.11133*alph2)
     e1=(-0.50003 - 0.18329*alph2 + 0.23835*alph2^2)/(1 + 1.6756*alph2)
     e2=(-0.50001 + 1.1414*alph2 + 0.0083589*alph2^2)/(1 + 11.168*alph2)
     e3=(-0.70003 - 0.73587*alph2 + 0.51509*alph2^2)/(1 + 4.7665*alph2)
     ua2 = log(1 + ama*AB2)    #below Eq.7a
     lam2 = (e0 + e1*ua2 + e2*ua2^2)/(1 + e3*ua2)
     
     #Aeroso transmittance   #below Eq.7a
     ta2=AB2*lam2^(-alph2)
     TA2=exp(-ama*ta2)
     TAS2=exp(-ama*0.84*ta2)   #w2=0.84 recommended by author
     
     #forward scattering fractions for Rayleigh extinction
     BR2=0.5  # multi scatter negibile in Band 2
     #the aerosol forward scatterance factor
     Ba=1 - exp(-0.6931-1.8326*cos(sza))

     #Aerosol scattering correction
     h0=(3.4352 + 0.65267*ama + 0.00034328*ama^2)/(1 + 0.034388*ama^1.5)
     h1=(1.231 - 1.63853*ama + 0.20667*ama^2)/(1 + 0.1451*ama^1.5)
     h2=(0.8889 - 0.55063*ama + 0.50152*ama^2)/(1 + 0.14865*ama^1.5)
     F2=(h0 + h1*ta2)/(1 + h2*ta2)

     #sky albedo     Appendix
     rs2=(0.010191 + 0.00085547*alph2 + AB2*(0.14618 + 0.062758*alph2)/(1 - 0.19402*alph2))/(1 + AB2*(0.58101 + 0.17426*alph2)/(1 - 0.17586*alph2))
     #ground albedo
     rg=albedo 
     
     #irradiance BAND1
     E0n1=Eext*0.46512
     #direct beam irradiance
     Ebn1=E0n1*TR1*Tg1*To1*Tn1*Tw1*TA1  
     #the incident diffuse irradiance on a perfectly absorbing ground 
     Edp1=E0n1*cos(sza)*To1*Tg1*Tn1166*Tw1166*(BR1*(1-TR1)*TA1^0.25 + Ba*F1*TR1*(1-TAS1^0.25))
     #multiple reflections between the ground and the atmosphere    
     Edd1=rg*rs1*(Ebn1*cos(sza)+Edp1)/(1-rg*rs1)   
     
     #irradiance BAND2
     E0n2 = Eext*0.51951
     #direct beam irradiance
     Ebn2 = E0n2*TR2*Tg2*To2*Tn2*Tw2*TA2
	 #the incident diffuse irradiance on a perfectly absorbing ground 
     Edp2 = E0n2*cos(sza)*To2*Tg2*Tn2166*Tw2166*(BR2*(1-TR2)*TA2^0.25 + Ba*F2*TR2*(1-TAS2^0.25))
     #multiple reflections between the ground and the atmosphere
     Edd2 = rg*rs2*(Ebn2*cos(sza)+Edp2)/(1-rg*rs2)

     #TOTALS BAND1+BAND2
     #direct norm irradiance
     Ebnrest2=(Ebn1 + Ebn2) 
     #diffuse horizontal irradiance
     Edhrest2=Edp1 + Edd1 + Edp2 + Edd2
     #global horizontal irradiance
     Eghrest2 = Ebnrest2* cos(sza) + Edhrest2

	###Quality control
    lower=0
    Ebnrest2[Ebnrest2<lower]=0
    Edhrest2[Edhrest2<lower]=0
    Eghrest2[Eghrest2<lower]=0
	return(list(Ebnrest2, Edhrest2, Eghrest2))
	}


