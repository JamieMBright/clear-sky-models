#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 35-40 Heliosat-2 2002
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Lefevre, M., Albuisson, M., & Wald, L. (2002). Joint report on interpolation scheme ‘meteosat’and database ‘climatology’i (meteosat). SoDa Deliverable D3-8 and D5-1-4. Internal document.

###Inputs:
#  Esc=1367   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  altitude   [m]         (location's altitude)
#  TL2        [double]    (Linke Turbidity at air mass = 2)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceHeliosat2<-function(altitude){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
	
	#air mass
	#correct alpha(solar altitude angle, alpha=pi/2-sza) for refraction
	alpha=pi/2-sza
	alphatrue=alpha+0.061359*(0.1594+1.123*alpha+0.065656*alpha^2)/(1+28.9344*alpha+277.3971*alpha^2)   
	am=1/(sin(alphatrue)+0.50572*(alphatrue/pi*180+6.07995)^-1.6364)
	#corr_rot is the correction of the integral Rayleigh optical thickness due to the elevation of the site 	
	corr_rot1=1
	corr_rot075=1.248174-0.011997*am+0.00037*am^2
	corr_rot05=1.68219-0.03059*am+0.00089*am^2
	
	#piecewise linear interpolation
	ntotal=length(corr_rot05)
	corr_rot=corr_rot05   #define the length of corr_rot
	xzuobiao=c(0.5,0.75,1)
	pchup0=exp(-altitude/8434.5)     #all sites satisfy p/p0>=0.5
	for (i in 1:ntotal){
		if (is.na(corr_rot05[i])&is.na(corr_rot075[i])) next
		corr_rot[i]=approx(xzuobiao,c(corr_rot05[i],corr_rot075[i],1),xout=pchup0)$y
	}
	
	#Rayleigh Optical Thickness (rot)    
    rot = (corr_rot^-1)*(6.625928+1.92969*am-0.170073*am^2+0.011517*am^3-0.000285*am^4)^-1
    rot[which(am>20)]=(10.4 + 0.718*am[which(am>20)]*exp(-altitude/8434.5))^-1
	
	L00=alpha    #define the length of L00 to be the same as the length of alpha
	L00[which(alpha>pi/6)]=-1.7349e-2
	L00[which((alpha>pi/12)&(alpha<=pi/6))]=-8.2193e-3
	L00[which(alpha<=pi/12)]=-1.1656e-3
    L01=alpha    #define the length
	L01[which(alpha>pi/6)]=-5.8985e-3
	L01[which((alpha>pi/12)&(alpha<=pi/6))]=4.5643e-4
	L01[which(alpha<=pi/12)]=1.8408e-4
	L02=alpha   #define the length
	L02[which(alpha>pi/6)]=6.8868e-4
	L02[which((alpha>pi/12)&(alpha<=pi/6))]=6.7916e-5
	L02[which(alpha<=pi/12)]=-4.8754e-7
   
    L10=alpha    #define the length
	L10[which(alpha>pi/6)]=1.0258
	L10[which((alpha>pi/12)&(alpha<=pi/6))]=8.9233e-1
	L10[which(alpha<=pi/12)]=7.4095e-1
    L11=alpha    #define the length
	L11[which(alpha>pi/6)]=-1.2196e-1
	L11[which((alpha>pi/12)&(alpha<=pi/6))]=-1.9991e-1
	L11[which(alpha<=pi/12)]=-2.2427e-1
	L12=alpha    #define the length
	L12[which(alpha>pi/6)]=1.9299e-3
	L12[which((alpha>pi/12)&(alpha<=pi/6))]=9.9741e-3
	L12[which(alpha<=pi/12)]=1.5314e-2
	
	L20=alpha    #define the length
	L20[which(alpha>pi/6)]=-7.2178e-3
	L20[which((alpha>pi/12)&(alpha<=pi/6))]=2.5428e-1
	L20[which(alpha<=pi/12)]=3.4959e-1
    L21=alpha    #define the length
	L21[which(alpha>pi/6)]=1.3086e-1
	L21[which((alpha>pi/12)&(alpha<=pi/6))]=2.6140e-1
	L21[which(alpha<=pi/12)]=7.2313e-1
	L22=alpha    #define the length
	L22[which(alpha>pi/6)]=-2.8405e-3
	L22[which((alpha>pi/12)&(alpha<=pi/6))]=-1.7020e-2
	L22[which(alpha<=pi/12)]=-1.2305e-1
	L23=alpha    #define the length
	L23[which(alpha>pi/6)]=0
	L23[which((alpha>pi/12)&(alpha<=pi/6))]=0
	L23[which(alpha<=pi/12)]=5.9194e-3
		
	#Trb=exp(-0.8662*TL2*exp(altitude/8434.5)*rot)
	
	#c0=L00+L01*TL2*pchup0+L02*(TL2*pchup0)^2
	#c1=L10+L11*TL2*pchup0+L12*(TL2*pchup0)^2
	#c2=L20+L21*TL2*pchup0+L22*(TL2*pchup0)^2+L23*(TL2*pchup0)^3
	#Fb=c0+c1*sin(alpha)+c2*(sin(alpha))^2
	
	#direct horizontal irradiance
	#Eb=Eext*pmax(Trb*Fb,0)
	
	#the diffuse transmission function at zenith
	#Trd=(-1.5843e-2)+(3.0543e-2)*TL2*pchup0+(3.797e-4)*(pchup0*TL2)^2
	#the diffuse angular function
    #a0 = (2.6463e-1)-(6.1581e-2)*TL2*pchup0+(3.1408e-3)*(pchup0*TL2)^2
    #a1 = 2.0402+(1.8945e-2)*TL2*pchup0-(1.1161e-2)*(pchup0*TL2)^2
    #a2 = -1.3025+(3.9231e-2)*TL2*pchup0+(8.5079e-3)*(pchup0*TL2)^2
    #a0[which((a0*Trd<(2e-3))==T)] = (2e-3)/Trd[which((a0*Trd<(2E-3))==T)]
    #FD = a0+a1*sin(alpha)+a2*(sin(alpha))^2
    #diffuse horizontal irradiance
    #Ed = Eext*pmax(Trd*FD,0)
    #global horizontal irradiance
    #Eghheliosat2=Eb+Ed
    
    TL2=list(TL2R,TL2D,TL2I,TL2Gu,TL2M,TL2Gr)
    EghHeliosat2_name=c('EghHeliosat2_R','EghHeliosat2_D','EghHeliosat2_I','EghHeliosat2_Gu','EghHeliosat2_M','EghHeliosat2_Gr')
    EbnHeliosat2_name=c('EbnHeliosat2_R','EbnHeliosat2_D','EbnHeliosat2_I','EbnHeliosat2_Gu','EbnHeliosat2_M','EbnHeliosat2_Gr')
    EdhHeliosat2_name=c('EdhHeliosat2_R','EdhHeliosat2_D','EdhHeliosat2_I','EdhHeliosat2_Gu','EdhHeliosat2_M','EdhHeliosat2_Gr')
    EghHeliosat2=vector('list',6)
    EdhHeliosat2=vector('list',6)
    EbnHeliosat2=vector('list',6)
    lower=0
    
    for (i in 1:6){   
    	#direct beam irradiance
    	Trb=exp(-0.8662*TL2[[i]]*exp(altitude/8434.5)*rot)
    	c0=L00+L01*TL2[[i]]*pchup0+L02*(TL2[[i]]*pchup0)^2
    	c1=L10+L11*TL2[[i]]*pchup0+L12*(TL2[[i]]*pchup0)^2
    	c2=L20+L21*TL2[[i]]*pchup0+L22*(TL2[[i]]*pchup0)^2+L23*(TL2[[i]]*pchup0)^3
    	Fb=c0+c1*sin(alpha)+c2*(sin(alpha))^2
    	EbnHeliosat2[[i]]=Eext*pmax(Trb*Fb,0)/sin(alpha)
    	EbnHeliosat2[[i]][EbnHeliosat2[[i]]<lower]=0       #QC
    	assign(EbnHeliosat2_name[i], EbnHeliosat2[[i]])
    	
    	#diffuse horizontal irradiance
    	TRD = (-1.5843e-2)+(3.0543e-2)*TL2[[i]]+(3.797e-4)*(TL2[[i]])^2
    	a01 = (2.6463e-1)-(6.1581e-2)*TL2[[i]]+(3.1408e-3)*(TL2[[i]])^2
    	a11 = 2.0402+(1.8945e-2)*TL2[[i]]-(1.1161e-2)*(TL2[[i]])^2
    	a21 = -1.3025+(3.9231e-2)*TL2[[i]]+(8.5079e-3)*(TL2[[i]])^2
    	xiabiao=which((a01*TRD<(2e-3))==T)
    	a01[xiabiao] = (2e-3)/TRD[xiabiao]
    	FD = a01+a11*sin(alpha)+a21*(sin(alpha))^2
    	EdhHeliosat2[[i]] = Eext*pmax(0,TRD*FD)
    	EdhHeliosat2[[i]][EdhHeliosat2[[i]]<lower]=0      #QC
    	assign(EdhHeliosat2_name[i], EdhHeliosat2[[i]])
    	
    	#global horizontal irradiance
    	EghHeliosat2[[i]]=EbnHeliosat2[[i]]*cos(sza)+EdhHeliosat2[[i]]
    	EghHeliosat2[[i]][EghHeliosat2[[i]]<lower]=0     #QC
    	assign(EghHeliosat2_name[i], EghHeliosat2[[i]])
    	}
 return(list(EghHeliosat2_R,EghHeliosat2_D,EghHeliosat2_I,EghHeliosat2_Gu,EghHeliosat2_M,EghHeliosat2_Gr,EbnHeliosat2_R,EbnHeliosat2_D,EbnHeliosat2_I,EbnHeliosat2_Gu,EbnHeliosat2_M,EbnHeliosat2_Gr, EdhHeliosat2_R,EdhHeliosat2_D,EdhHeliosat2_I,EdhHeliosat2_Gu,EdhHeliosat2_M,EdhHeliosat2_Gr))
 }

