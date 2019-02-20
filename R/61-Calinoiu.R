#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 61-Calinoiu 2018
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Calinoiu, D., Stefu, N., Boata, R., Blaga, R., Pop, N., Paulescu, E., ... & Paulescu, M. (2018). Parametric modeling: A simple and versatile route to solar irradiance. Energy Conversion and Management, 164, 175-187.

###Inputs:
#  Esc=1366.1   [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
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
IrradianceCalinoiu<-function(){	
	#extraterrestrial irradiance
	Esc=1366.1
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
	
	#Beam normal irradiance
	#ozone
	mo3=1/(cos(sza)+268.45*((sza/pi*180)^0.5)*(115.42-sza/pi*180)^(-3.2922) )
	xo3=mo3*ozone
	To3= (1+8.5951*xo3+0.2179*xo3^2) / (1+8.75308*xo3+0.45*xo3^2-0.0004*xo3^3)
	
	#nitrogen dioxide
	mno2=1/(cos(sza)+602.3*((sza/pi*180)^0.5)*(117.96-sza/pi*180)^(-3.4536) )
	xno2=mno2*NO2
	Tno2=0.742+0.258*exp(-xno2/0.09)
	
	#water vapour
	mw=1/(cos(sza)+0.031141*((sza/pi*180)^0.1)*(92.471-sza/pi*180)^(-1.3814) )
	xw=mw*wv
	Tw=(1+1.4*xw^0.5+0.053*xw)*(1-0.013*xw)/(1+1.626414*xw^0.5+0.10816267*xw)
	
	#trace gas absorption
	mg=1/(cos(sza)+0.45665*((sza/pi*180)^0.07)*(96.4836-sza/pi*180)^(-1.6970) )
	xg=mg
	Tg=(1+0.19558*xg^0.5)/(1+0.215582*xg^0.5+0.0005*xg)
	
	#Rayleigh scattering
	mR=1/(cos(sza)+0.45665*((sza/pi*180)^0.07)*(96.4836-sza/pi*180)^(-1.6970) )
	xR=mR
	TR=(1+0.1564*xR+0.0001*xR^2)/(1+0.26038*xR+0.00697*xR^2)
	
	#aerosol extinction
	ma=1/(cos(sza)+0.031141*((sza/pi*180)^0.1)*(92.471-sza/pi*180)^(-1.3814) )
	xa=ma*ang_beta
	Ta=(1-0.046*xa)/(1+1.73849*xa+0.79081*xa^2)
	
	taub=To3*Tno2*Tw*Tg*TR*Ta
	
	#direct beam irradiance
	EbnCalinoiu=Eext*taub
	
	#Diffuse irradiance	
	#aerosol scattering
	g=0.5
	f1g05=1.459*log(1-g)+0.1595*(log(1-g))^2+0.4129*(log(1-g))^3
	f2g05=0.0783*log(1-g)-0.3824*(log(1-g))^2+0.5874*(log(1-g))^3
	gamma_a05=0.92*(1-0.5*exp(f1g05*cos(sza)+f2g05*(cos(sza))^2))
	g=0.75
	f1g075=1.459*log(1-g)+0.1595*(log(1-g))^2+0.4129*(log(1-g))^3
	f2g075=0.0783*log(1-g)-0.3824*(log(1-g))^2+0.5874*(log(1-g))^3
	gamma_a075=0.92*(1-0.5*exp(f1g075*cos(sza)+f2g075*(cos(sza))^2))
	gamma_a=0.5*(gamma_a05+gamma_a075)
	
	omega_hat=0.95
	Tas=(0.96865-1.57909*log(omega_hat)-0.07119*xa+2.78513*(log(omega_hat))^2+0.00157*xa^2+0.07939*xa*log(omega_hat))/(1-1.22067*log(omega_hat)+1.83599*xa+3.20558*(log(omega_hat))^2+0.55453*xa^2-1.73402*xa*log(omega_hat))
	Taa=Ta/Tas*(0.14578+0.69109/omega_hat)
	
	mo3[mo3>=17]=16.999999999
	mo3[mo3<=0.001]=0.001
    xo3[xo3>=0.45]=0.449999999
	xo3[xo3<=0.15]=0.150000001
	xo3_1=-1+2/0.3*(xo3-0.15)
	mo3_1=-1+2/17*(mo3-0.01)
	xo3_1[xo3_1==1]=0.9999999
	mo3_1[mo3_1==1]=0.9999999
	cij=matrix(0,nrow=10,ncol=10)
	cij[1,]=c(0.907, -0.031, 1.026E-3, -1.106E-4, 1.703E-5, -4.738E-6, 4.409E-6, 1.145E-6, -1.49E-6, 1.748E-8)
	cij[2,]=c(-0.044, -0.014, 8.96E-4, -8.85E-5, 5.172E-6, -3.739E-6, -7.45E-7, -3.646E-6, 1.707E-6, 7.464E-8)
	cij[3,]=c(0.023, 7.475E-3, -3.653e-4, 3.38E-5, -1.041e-6, -4.724e-6, -8.628e-7, 2.768e-7, 2.305e-6, 0)
	cij[4,]=c(-0.012, -4.014E-3, 1.514E-4, -1.806E-5, 2.288E-6, -2.995E-6, -3.292E-6, -1.753E-6, 0, 0)
	cij[5,]=c(6.22E-3, 2.046E-3, -5.799E-5, 8.207E-6, -1.759E-6, -1.807E-6, 8.908E-7, 0, 0, 0)
	cij[6,]=c(3.158E-3, -9.618E-4, 2.104E-5, -4.049E-6, 3.605E-6, -1.828E-6, 0, 0, 0, 0)
	cij[7,]=c(1.576E-3, 4.2E-4, -1.083E-5, 7.229E-6, -3.751E-9, 0, 0, 0, 0, 0)
	cij[8,]=c(-8.093E-4, -1.768E-4, 6.057E-6, -1.74E-6, 0, 0, 0, 0, 0, 0)
	cij[9,]=c(4.428E-4, 6.802E-5, -4.685E-6, 0, 0, 0, 0, 0, 0, 0)
	cij[10,]=c(-2.212E-4, -3.167E-5, 0, 0, 0, 0, 0, 0, 0, 0)
	Td_o3=0
	for (i in 0:9){
		for (j in 0:(9-i)){
			Xjxo3=cos(j*acos(xo3_1))
			Ximo3=cos(i*acos(mo3_1))
			Td_o3=Td_o3+cij[(i+1),(j+1)]*Xjxo3*Ximo3
		}
	}
	
	tauda=gamma_a*Td_o3*Tno2*Tw*Tg*TR*Taa*(1-Tas)
	
	#Rayleigh
	taudR=0.5*Td_o3*Tno2*Tw*Tg*Taa*(1-TR)
	
	#diffuse horizontal irradiance
	EdhCalinoiu=Eext*cos(sza)*(taudR+tauda)
	
	#Global irradiance
	EghCalinoiu=EbnCalinoiu*cos(sza)+EdhCalinoiu
	
	###Quality control
    lower=0
    EbnCalinoiu[EbnCalinoiu<lower]=0
    EdhCalinoiu[EdhCalinoiu<lower]=0
    EghCalinoiu[EghCalinoiu<lower]=0
	return(list(EbnCalinoiu, EdhCalinoiu, EghCalinoiu))
	}

