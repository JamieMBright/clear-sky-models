#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 54-CEM 1978
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Atwater, M. A., & Ball, J. T. (1978). A numerical solar radiation model based on standard meteorological observations. Solar Energy, 21(3), 163-170.
#Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

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
IrradianceCEM<-function(){	
	#extraterrestrial irradiance
	Esc=1353
    #Eext=Esc*(1+0.033*cos(2*pi*Dayth/365))
    totaldayofyear=366-ceiling(Year/4-trunc(Year/4))  
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    	
    am=35/(1224*(cos(sza))^2+1)^0.5 
    amp=am*press/1013.25
    p=press/1013.25*101.325
    X1=(am*(949*p*1e-5+0.051))^0.5 
    Trb=1.041-0.16*X1
    Trg=1.021-0.0824*X1 
    Tg=1
    To=1
    absw=0.077*(am*wv)^0.3
    Trw=1-absw
    Tra=exp(-amp*broadbandaod)
    
    #direct normal irradiance
    EbnCEM=Eext*Tra*(Trb*Tg-absw)
    #global horizontal irradiance
    EghCEM=Eext*cos(sza)*Tra*(Trg*Tg-absw)/(1-albedo*0.0685) 
    #diffuse horizontal irradiance
    EdhCEM=EghCEM-EbnCEM*cos(sza) 
     
    ###Quality control
    lower=0
    EbnCEM[EbnCEM<lower]=0
    EdhCEM[EdhCEM<lower]=0
    EghCEM[EghCEM<lower]=0
	return(list(EbnCEM,EdhCEM,EghCEM))
	}

