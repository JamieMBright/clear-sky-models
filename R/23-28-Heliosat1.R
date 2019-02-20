#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 23-28 Heliosat-1 1996
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Hammer, A., Heinemann, D., Hoyer, C., Kuhlemann, R., Lorenz, E., Müller, R., & Beyer, H. G. (2003). Solar energy assessment using remote sensing technologies. Remote Sensing of Environment, 86(3), 423-432.
#Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

###Inputs:
#  Esc=1367   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  press      [mb]        (local barometric)
#  TL2        [double]    (Linke Turbidity at air mass = 2)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceHeliosat1<-function(){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))

	#absolute air mass
	AM=1/(cos(sza)+0.50572*((6.07995+(pi/2-sza)*180/pi)^(-1.6364)))	
	ama=AM*press/1013.25
	
	#the Rayleigh optical thickness
	theta=1/(6.62960+1.7513*ama-0.1202*ama^2+0.0065*ama^3-0.00013*ama^4)
	theta[which(ama>20)]=(10.4 + 0.718*ama[which(ama>20)])^-1
	
	#direct beam irradiance
	#EbnHeliosat1=Eext*exp(-ama*theta*TL2*0.8662)
	#diffuse horizontal irradiance
	#EdhHeliosat1=Eext*(0.0065+(0.0646*TL2-0.045)*cos(sza)-(0.0327*TL2-0.014)*(cos(sza))^2)
	#global horizontal irradiance
    #EghHeliosat1=EbnHeliosat1*cos(sza)+EdhHeliosat1

    TL2=list(TL2R,TL2D,TL2I,TL2Gu,TL2M,TL2Gr)
    EghHeliosat1_name=c('EghHeliosat1_R','EghHeliosat1_D','EghHeliosat1_I','EghHeliosat1_Gu','EghHeliosat1_M','EghHeliosat1_Gr')
    EbnHeliosat1_name=c('EbnHeliosat1_R','EbnHeliosat1_D','EbnHeliosat1_I','EbnHeliosat1_Gu','EbnHeliosat1_M','EbnHeliosat1_Gr')
    EdhHeliosat1_name=c('EdhHeliosat1_R','EdhHeliosat1_D','EdhHeliosat1_I','EdhHeliosat1_Gu','EdhHeliosat1_M','EdhHeliosat1_Gr')
    EghHeliosat1=vector('list',6)
    EdhHeliosat1=vector('list',6)
    EbnHeliosat1=vector('list',6)
    lower=0
    
    #calculate Ebn,Edh,Egh and quality control
    for (i in 1:6){
    	#direct beam irradiance
    	EbnHeliosat1[[i]]=pmax(0,(Eext*exp(-ama*theta*TL2[[i]]*0.8662)))
    	EbnHeliosat1[[i]][EbnHeliosat1[[i]]<lower]=0    #QC
    	assign(EbnHeliosat1_name[i],EbnHeliosat1[[i]])
    	
    	#diffuse horizontal irradiance
    	EdhHeliosat1[[i]]=pmax(0,Eext*(0.0065+(0.0646*TL2[[i]]-0.045)*cos(sza)-(0.0327*TL2[[i]]-0.014)*(cos(sza))^2))
    	EdhHeliosat1[[i]][EdhHeliosat1[[i]]<lower]=0     #QC
    	assign(EdhHeliosat1_name[i],EdhHeliosat1[[i]])
    	
    	#global horizontal irradiance
    	EghHeliosat1[[i]]=EbnHeliosat1[[i]]*cos(sza)+EdhHeliosat1[[i]]
    	EghHeliosat1[[i]][EghHeliosat1[[i]]<lower]=0      #QC
    	assign(EghHeliosat1_name[i],EghHeliosat1[[i]])
    }
	return(list(EghHeliosat1_R,EghHeliosat1_D,EghHeliosat1_I,EghHeliosat1_Gu,EghHeliosat1_M,EghHeliosat1_Gr,EbnHeliosat1_R,EbnHeliosat1_D,EbnHeliosat1_I,EbnHeliosat1_Gu,EbnHeliosat1_M,EbnHeliosat1_Gr, EdhHeliosat1_R,EdhHeliosat1_D,EdhHeliosat1_I,EdhHeliosat1_Gu,EdhHeliosat1_M,EdhHeliosat1_Gr))
	}
	
	