#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 41-46 Ineichen & Perez 2002
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Ineichen, P., & Perez, R. (2002). A new airmass independent formulation for the Linke turbidity coefficient. Solar Energy, 73(3), 151-157.

###Inputs:
#  Esc=1367   [Wm-2]      (Solar constant)
#  sza        [radians]   (zenith_angle) 
#  altitude   [m]         (location's altitude)
#  TL2        [double]    (Linke Turbidity at air mass = 2)

###Outputs:
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceIp<-function(altitude){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))  
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #air mass
    AM=1/(cos(sza)+0.50572*((6.07995+(pi/2-sza)*180/pi)^(-1.6364)))
    
	#coefficients based on altitude
	fh1 = exp(-altitude/8000)
    fh2 = exp(-altitude/1250)
    cg1 = (5.09e-5)*(altitude)+0.868
    cg2 = (3.92e-5)*(altitude)+0.0387
    
    #global horizontal irradiance & QC
    #Eghip = cg1*Eext*cos(sza)*exp(-cg2*AM*(fh1+fh2*(TL2-1)))
    
    TL2=list(TL2R,TL2D,TL2I,TL2Gu,TL2M,TL2Gr)
    Eghip_name=c('Eghip_R','Eghip_D','Eghip_I','Eghip_Gu','Eghip_M','Eghip_Gr')
    Eghip=vector('list',6)
    lower=0
    for (i in 1:6){    
    	Eghip[[i]]=pmax(0,(cg1*Eext*cos(sza)*exp(-cg2*AM*(fh1+fh2*(TL2[[i]]-1)))))
    	Eghip[[i]][Eghip[[i]]<lower]=0   #QC
    	assign(Eghip_name[i],Eghip[[i]])
    }
	return(list(Eghip_R,Eghip_D,Eghip_I,Eghip_Gu,Eghip_M,Eghip_Gr))
	}


