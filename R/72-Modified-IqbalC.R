#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 71-Modified Iqbal-C 2012
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Bird, R. E., & Hulstrom, R. L. (1981). Simplified clear sky model for direct and diffuse insolation on horizontal surfaces (No. SERI/TR-642-761). Solar Energy Research Inst., Golden, CO (USA).
#Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

###Inputs:
#  Esc=1367     [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  altitude     [m]                (location's altitude)
#  albedo       [double]           (surface/ground albedo)
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
IrradianceMIqbalc<-function(altitude){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))

    #air mass
    am=(cos(sza)+0.15*(3.885+(pi/2-sza)/pi*180)^-1.253)^-1
    ama=am*exp(-0.0001184*altitude)     #1st difference by Gueymard2012
   
    #the transmittance by Rayleigh scattering
    TR=exp(-0.0903*(ama)^0.84*(1+ama-ama^1.01)) 

     #the transmittance by ozone
     U3=ozone*am       
     alpho=(0.1611*U3*(1+139.48*U3)^-0.3035 - 0.002715*U3*(1.0+0.044*U3+0.0003*U3^2)^-1)
     TO=1-alpho

     #the transmittance by uniformly mixed gases
     TG=exp(-0.0127*ama^0.26)

     #the transmittance by water vapor
     PW=wv*am   #2nd difference
     alphw=2.4959*PW*((1+79.034*PW)^0.6828+6.385*PW)^-1
     TW=1-alphw

     #the transmittance by Aerosol
     TA=(0.12445*ang_alpha-0.0162)+(1.003-0.125*ang_alpha)*exp(-ama*ang_beta*(1.089*ang_alpha+0.5123))    #3rd difference
     
     #direct beam irradiance
     EbnMiqbalc=0.9751*Eext*TR*TO*TG*TW*TA

     #the transmittance of direct radiation due to aerosol absorptance
     w0=0.9 # according to Bird Hulstrom
     TAA=1-(1-w0)*(1-ama+ama^1.06)*(1- TA)
     #the Rayleigh-scattered diffuse irradiance after the first pass through the atmosphere
     Edr = 0.79*Eext*cos(sza)*TO*TG*TW*TAA*0.5*(1-TR)/(1-ama+ama^1.02)
     
     #the aerosol scattered diffuse irradiance after the first pass through the atmosphere
     TAS=TA/TAA
     Fc=0.84 # according to Iqbal
     Eda=0.79*Eext*cos(sza)*TO*TG*TW*TAA*Fc*(1-TAS)/(1-ama+ama^1.02)

     #global horizontal irradiance
     rg=albedo   #ground albedo
     ra=0.0685+(1-Fc)*(1-TAS)
     EghMiqbalc=(EbnMiqbalc*cos(sza)+Edr+Eda)/(1-rg*ra)
     #diffuse horizontal irradiance
     EdhMiqbalc=EghMiqbalc-EbnMiqbalc*cos(sza)

	###Quality control
    lower=0
    EbnMiqbalc[EbnMiqbalc<lower]=0
    EdhMiqbalc[EdhMiqbalc<lower]=0
    EghMiqbalc[EghMiqbalc<lower]=0
	return(list(EbnMiqbalc, EdhMiqbalc, EghMiqbalc))
	}


