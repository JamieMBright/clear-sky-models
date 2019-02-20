#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 63-MAC2 1982
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Davies, J. A., & McKay, D. C. (1982). Estimating solar irradiance and components. Solar Energy, 29(1), 55-64.

###Inputs:
#  Esc=1353     [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  press        [mb]               (local barometric)
#  albedo       [double]           (surface/ground albedo)
#  ang_alpha    [dimensionless]    (Angstrom_exponent, also known as alpha)
#  ang_beta     [dimensionless]    (Angstrom turbidity coefficient)
#  wv           [atm.cm]           (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceMac2<-function(){
	#extraterrestrial irradiance
	Esc=1353
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
	
	#Air Mass
    amm=35/((1224*(cos(sza))^2+1)^0.5)
    amm[amm<0.0]=NA  

    #Ozone Transmittance 
    ozone=0.35   #Davies and Mckay 1982 set ozone a fixed value of 3.5mm    
    XO=amm*(ozone*10)    #Davies and Mckay 1982 ozone unit is mm, here in the code unit is cm
    AO=((0.1082*XO)/((1+13.86*XO)^0.805)) + ((0.00658*XO)/(1+(10.36*XO)^3)) + (0.002118*XO/(1+0.0042*XO+3.23e-6*XO^2))
    TO=1-AO

    #Rayleigh Transmittance  (linear interpolation based on  Table 2 in Davies and Mckay 1982 )
    TR=amm
    amms=c(0.5,1,1.2,1.4,1.6,1.8,2.0,3.0,3.5,4.0,4.5,5.0,5.5,6,10,30)
    TRs=c(.9385,.8973,.8830,.8696,.8572,.8455,.8344,.7872,.7673,.7493,.7328,.7177,.7037,.6907,.6108,.4364)
    TR=approx(amms,TRs,xout=amm)$y
    #TR = 0.9768-0.0874*amm + 0.010607552*amm^2 - (8.46205E-4)*amm^3	+ (3.57246E-5)*amm^4 - (6.0176E-7)*amm^5  (Davies et al. 1988,P21)
    
    #Aerosols Transmittance, borrowed from BH 1981
    tA=ang_beta*(0.38^-ang_alpha*0.2758+0.5^-ang_alpha*0.35)
    TA=exp(-tA*amm) 

    #Water Vapor Transmittance
    XW=amm*wv*10*(press/1013.25)^0.75     #we do not use dewpoint tem to calculate wv, for we have wv already
    AW=0.29*XW/(((1+ 14.15*XW)^0.635)+0.5925*XW)

    #Forward Scatter
    szajiao=sza*180/pi
    szajiaos=c(0,25.8,36.9,45.6,53.1,60.0,66.4,72.5,78.5,90)
    fs=c(.92,.91,.89,.86,.83,.78,.71,.67,.60,.60)
    f=szajiao
    f=approx(szajiaos,fs,xout=szajiao)$y
    
    #direct beam irradiance
    Ebnmac=Eext*(TO*TR-AW)*TA
    
    #diffuse components from Rayleigh scatter
    DR=Eext*cos(sza)*TO*(1-TR)/2
    #diffuse components from scattering by aerosol
    DA=Eext*cos(sza)*(TO*TR-AW)*(1-TA)*0.75*f     # w0 = 0.75 according to Table5 in Davies and Mckay 1982
    #diffuse horizontal irradiance
    Taaa=0.95^1.66     #Taaa is TA determined at amm=1.66 ,   k=0.95 according to Table5 in Davies and Mckay 1982
    poub=0.0685+(1-Taaa)*0.75*(1-0.83)   #f' is f determined at amm=1.66, f' equals  0.83, estimate theta when amm=1.66, theta near 53 degree
    Edhmac=poub*albedo*(Ebnmac*cos(sza)+DR+DA)/(1-poub*albedo)+DR+DA
    
    #global horizontal irradiance
    Eghmac=(Ebnmac*cos(sza)+DR+DA)/(1-poub*albedo)

    ###Quality control
    lower=0
    Ebnmac[Ebnmac<lower]=0
    Edhmac[Edhmac<lower]=0
    Eghmac[Eghmac<lower]=0
	return(list(Ebnmac, Edhmac, Eghmac))
	}
	
	
