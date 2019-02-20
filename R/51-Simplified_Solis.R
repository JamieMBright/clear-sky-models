#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 51-Simplified Soils 2008
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Ineichen, P. (2008). A broadband simplified version of the Solis clear sky model. Solar Energy, 82(8), 758-762.

###Inputs:
#  Esc=1367   [Wm-2]              (Solar constant)
#  sza        [radians]           (zenith_angle) 
#  press      [mb]                (local barometric)
#  AOD_700    [dimensionless]     (aerosol optical depth at 700 nm)
#  wv         [atm.cm]            (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceSolis<-function(){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #Ineichen 2008 limit aod700 between 0~0.45    
    aod700=pmin(0.45,aod700)
    
	#Modified solar constant(Eext)
    Eext0 = 1.08*wv^0.0051
    Eext1 = 0.97*wv^0.032
    Eext2 = 0.12*wv^0.56
    Eexte = Eext*(Eext2*aod700^2 + Eext1*aod700 + Eext0 + 0.071*log(press/1013.25))
 
    #direct total optical depth
    tb1 = 1.82 + 0.056*log(wv) + 0.0071*(log(wv))^2
    tb0 = 0.33 + 0.045*log(wv) + 0.0096*(log(wv))^2
    tbp = 0.0089*wv + 0.13
    tb = tb1*aod700 + tb0 + tbp*log(press/1013.25)
    b1= 0.00925*aod700^2+0.0148*aod700 - 0.0172
    b0= -0.7565*aod700^2+0.5057*aod700 + 0.4557
    b = b1*log(wv) + b0
    
    #diffuse total optical depth
    #for aod700<0.05
    td4 = 86*wv - 13800
    td3 = -3.11*wv + 79.4
    td2 = -0.23*wv + 74.8
    td1 = 0.092*wv - 8.86
    td0 = 0.0042*wv + 3.12
    tdp = -0.83*(1+aod700)^-17.2
    #for aod700>=0.05
    td4[which(aod700>=0.05)] = -0.21*wv[which(aod700>=0.05)] + 11.6
    td3[which(aod700>=0.05)] = 0.27*wv[which(aod700>=0.05)] - 20.7
    td2[which(aod700>=0.05)] = -0.134*wv[which(aod700>=0.05)] + 15.5
    td1[which(aod700>=0.05)] = 0.0554*wv[which(aod700>=0.05)] - 5.71
    td0[which(aod700>=0.05)] = 0.0057*wv[which(aod700>=0.05)] + 2.94
    tdp[which(aod700>=0.05)] = -0.71*(1+aod700[which(aod700>=0.05)])^-15
    td = td4*aod700^4 + td3*aod700^3 + td2*aod700^2 + td1*aod700 +td0 + tdp*log(press/1013.25)
    dp = 1/(18+152*aod700)
    d = -0.337*aod700^2 + 0.63*aod700 + 0.116 + dp*log(press/1013.25)
    
    #global total optical depth
    tg1 = 1.24 + 0.047*log(wv) + 0.0061*(log(wv))^2
    tg0 = 0.27 + 0.043*log(wv) + 0.0090*(log(wv))^2
    tgp = 0.0079*wv + 0.1
    tg = tg1*aod700 + tg0 + tgp*log(press/1013.25)
    g = -0.0147*log(wv) - 0.3079*aod700^2 + 0.2846*aod700 + 0.3798
    
    #direct normal radiation
	Ebnsolis=Eexte*exp(-tb/((cos(sza))^b))
	
	#diffuse horizontal radiation
	Edhsolis=Eexte*exp(-td/((cos(sza))^d))
	
    #global horizontal irradiance
    Eghsolis=Eexte*exp(-tg/((cos(sza))^g))*cos(sza)
    
	###Quality control
    lower=0
    Ebnsolis[Ebnsolis<lower]=0
    Edhsolis[Edhsolis<lower]=0
    Eghsolis[Eghsolis<lower]=0
	return(list(Ebnsolis, Edhsolis, Eghsolis))
	}


