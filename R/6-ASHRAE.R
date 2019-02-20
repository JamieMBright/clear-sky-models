#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 06-ASHRAE 1985
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Handbook, A. F. (1985). American society of heating, refrigerating and air-conditioning engineers. Inc.: Atlanta, GA, USA.
#El Mghouchi, Y., Ajzoul, T., Taoukil, D., & El Bouardi, A. (2016). The most suitable prediction model of the solar intensity, on horizontal plane, at various weather conditions in a specified location in Morocco. Renewable and Sustainable Energy Reviews, 54, 84-98.

###Inputs:
#  sza  [radians]  (zenith_angle) 

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Month is the month number ranging from 1 to 12. 

###Codes:
IrradianceASHRAE<-function(){
	#coeff of A,B,C from Jan,Feb,....,Dec
	A=c(1230,1215,1186,1136,1104,1088,1085,1107,1152,1193,1221,1234)
	B=c(0.142,0.144,0.156,0.180,0.196,0.205,0.207,0.201,0.177,0.160,0.149,0.142)
	C=c(0.058,0.060,0.071,0.097,0.121,0.134,0.136,0.122,0.092,0.073,0.063,0.057)
	#direct normal irradiance
	#EbnASHRAE=A*exp(-B/sinalpha) 
	EbnASHRAE=A[1]*exp(-B[1]/(cos(sza)))   #Jan
	EbnASHRAE[which(Month==2)]=A[2]*exp(-B[2]/(cos(sza))[which(Month==2)])  #Feb
	EbnASHRAE[which(Month==3)]=A[3]*exp(-B[3]/(cos(sza))[which(Month==3)])  #Mar
	EbnASHRAE[which(Month==4)]=A[4]*exp(-B[4]/(cos(sza))[which(Month==4)])  #Apr
	EbnASHRAE[which(Month==5)]=A[5]*exp(-B[5]/(cos(sza))[which(Month==5)])  #May
	EbnASHRAE[which(Month==6)]=A[6]*exp(-B[6]/(cos(sza))[which(Month==6)])  #Jun
	EbnASHRAE[which(Month==7)]=A[7]*exp(-B[7]/(cos(sza))[which(Month==7)])  #Jul
	EbnASHRAE[which(Month==8)]=A[8]*exp(-B[8]/(cos(sza))[which(Month==8)])  #Aug
	EbnASHRAE[which(Month==9)]=A[9]*exp(-B[9]/(cos(sza))[which(Month==9)])  #Sep
	EbnASHRAE[which(Month==10)]=A[10]*exp(-B[10]/(cos(sza))[which(Month==10)])  #Oct
	EbnASHRAE[which(Month==11)]=A[11]*exp(-B[11]/(cos(sza))[which(Month==11)])  #Nov
	EbnASHRAE[which(Month==12)]=A[12]*exp(-B[12]/(cos(sza))[which(Month==12)])  #Dec
	
	#diffuse horizontal irradiance
	#EdhASHRAE=C*EbnASHRAE
	EdhASHRAE=C[1]*EbnASHRAE   #Jan
	EdhASHRAE[which(Month==2)]=C[2]*EbnASHRAE[which(Month==2)]  #Feb
	EdhASHRAE[which(Month==3)]=C[3]*EbnASHRAE[which(Month==3)]  #Mar
	EdhASHRAE[which(Month==4)]=C[4]*EbnASHRAE[which(Month==4)]  #Apr
	EdhASHRAE[which(Month==5)]=C[5]*EbnASHRAE[which(Month==5)]  #May
	EdhASHRAE[which(Month==6)]=C[6]*EbnASHRAE[which(Month==6)]  #Jun
	EdhASHRAE[which(Month==7)]=C[7]*EbnASHRAE[which(Month==7)]  #Jul
	EdhASHRAE[which(Month==8)]=C[8]*EbnASHRAE[which(Month==8)]  #Aug
	EdhASHRAE[which(Month==9)]=C[9]*EbnASHRAE[which(Month==9)]  #Sep
	EdhASHRAE[which(Month==10)]=C[10]*EbnASHRAE[which(Month==10)]  #Oct
	EdhASHRAE[which(Month==11)]=C[11]*EbnASHRAE[which(Month==11)]  #Nov
	EdhASHRAE[which(Month==12)]=C[12]*EbnASHRAE[which(Month==12)]  #Dec
	
	#global horizontal irradiance
	EghASHRAE=EbnASHRAE*cos(sza)+EdhASHRAE
	
    ###Quality control
    lower=0
    EbnASHRAE[EbnASHRAE<lower]=0
    EdhASHRAE[EdhASHRAE<lower]=0
    EghASHRAE[EghASHRAE<lower]=0
	return(list(EbnASHRAE, EdhASHRAE, EghASHRAE))
	}
	
	
