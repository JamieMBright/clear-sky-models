#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 29-34 ESRA 2000
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Rigollier, C., Bauer, O., & Wald, L. (2000). On the clear sky model of the ESRA—European Solar Radiation Atlas—with respect to the Heliosat method. Solar energy, 68(1), 33-48.

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
IrradianceEsra<-function(altitude){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))     
	B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))

	#Air mass
	#correct alpha(solar altitude angle, alpha=pi/2-sza) for refraction
	alpha=pi/2-sza
	alphatrue=alpha+0.061359*(0.1594+1.123*alpha+0.065656*alpha^2)/(1+28.9344*alpha+277.3971*alpha^2)   
	ame=exp(-altitude/8434.5)/(sin(alphatrue)+0.50572*(alphatrue/pi*180+6.07995)^-1.6364)
    #Rayleigh Optical Thickness (rot)
    rot = (6.6296+1.7513*ame-0.1202*(ame^2)+0.0065*ame^3-0.00013*ame^4)^-1
    rot[which(ame>20)]=(10.4 + 0.718*ame[which(ame>20)])^-1
    
    #direct beam irradiance
    #Ebnesra = Eext*exp(-0.8662*TL2*ame*rot)
   
    #the diffuse transmission function at zenith
    #TRD = (-1.5843e-2)+(3.0543e-2)*TL2+(3.797e-4)*TL2^2
    #the diffuse angular function
    #a0 = (2.6463e-1)-(6.1581e-2)*TL2+(3.1408e-3)*TL2^2
    #a1 = 2.0402+(1.8945e-2)*TL2-(1.1161e-2)*TL2^2
    #a2 = -1.3025+(3.9231e-2)*TL2+(8.5079e-3)*TL2^2
    #a0[which((a0*TRD<(2e-3))==T)] = (2e-3)/TRD[which((a0*TRD<(2E-3))==T)]
    #FD = a0+a1*sin(alpha)+a2*(sin(alpha))^2
    #diffuse horizontal irradiance
    #Edesra = Eext*TRD*FD    
    
    #global horizontal irradiance
    #Eghesra= Ebnesra*cos(sza)+Edesra
    
    TL2=list(TL2R,TL2D,TL2I,TL2Gu,TL2M,TL2Gr)
    EghEsra_name=c('EghEsra_R','EghEsra_D','EghEsra_I','EghEsra_Gu','EghEsra_M','EghEsra_Gr')
    EbnEsra_name=c('EbnEsra_R','EbnEsra_D','EbnEsra_I','EbnEsra_Gu','EbnEsra_M','EbnEsra_Gr')
    EdhEsra_name=c('EdhEsra_R','EdhEsra_D','EdhEsra_I','EdhEsra_Gu','EdhEsra_M','EdhEsra_Gr')
    EghEsra=vector('list',6)
    EdhEsra=vector('list',6)
    EbnEsra=vector('list',6)
    lower=0
    for (i in 1:6){    #get('') can change character to variable name
    	#direct beam irradiance
    	EbnEsra[[i]]=pmax(0,Eext*exp(-0.8662*TL2[[i]]*ame*rot))
    	EbnEsra[[i]][EbnEsra[[i]]<lower]=0       #QC
    	assign(EbnEsra_name[i], EbnEsra[[i]])
    	#diffuse horizontal irradiance
    	TRD = (-1.5843e-2)+(3.0543e-2)*TL2[[i]]+(3.797e-4)*(TL2[[i]])^2
    	a01 = (2.6463e-1)-(6.1581e-2)*TL2[[i]]+(3.1408e-3)*(TL2[[i]])^2
    	a11 = 2.0402+(1.8945e-2)*TL2[[i]]-(1.1161e-2)*(TL2[[i]])^2
    	a21 = -1.3025+(3.9231e-2)*TL2[[i]]+(8.5079e-3)*(TL2[[i]])^2
    	xiabiao=which((a01*TRD<(2e-3))==T)
    	a01[xiabiao] = (2e-3)/TRD[xiabiao]
    	FD = a01+a11*sin(alpha)+a21*(sin(alpha))^2
    	EdhEsra[[i]] = Eext*TRD*FD
    	EdhEsra[[i]][EdhEsra[[i]]<lower]=0      #QC
    	assign(EdhEsra_name[i], EdhEsra[[i]])	
    	#global horizontal irradiance
    	EghEsra[[i]]=EbnEsra[[i]]*cos(sza)+EdhEsra[[i]]
    	EghEsra[[i]][EghEsra[[i]]<lower]=0      #QC
    	assign(EghEsra_name[i], EghEsra[[i]])
    }
    return(list(EghEsra_R,EghEsra_D,EghEsra_I,EghEsra_Gu,EghEsra_M,EghEsra_Gr,EbnEsra_R,EbnEsra_D,EbnEsra_I,EbnEsra_Gu,EbnEsra_M,EbnEsra_Gr, EdhEsra_R,EdhEsra_D,EdhEsra_I,EdhEsra_Gu,EdhEsra_M,EdhEsra_Gr))
	}


