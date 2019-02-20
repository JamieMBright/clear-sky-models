#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
### 52-Advanced Soils 2018
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----#
 
###References:
#Ineichen, P. (2018). High turbidity solis clear sky model: development and validation. Remote Sensing, 10(3), 435.
#Author's code

###Inputs:
#  Esc=1367   [Wm-2]              (Solar constant)
#  sza        [radians]           (zenith_angle) 
#  press      [mb]                (local barometric)
#  AOD_550    [dimensionless]     (aerosol optical depth at 550 nm)
#  wv         [atm.cm]            (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceSolis18<-function(){
	#extraterrestrial irradiance
	Esc=1367
	totaldayofyear=366-ceiling(Year/4-trunc(Year/4))
    B=(Dayth-1)*2*pi/totaldayofyear
    Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
    
    #settings
    En0=Eext
    aod=aod550
    aodln=log(aod)
    aod2=aod*aod
    aod3=aod*aod2
    aod4=aod*aod3
    aod5=aod*aod4
    pw=wv
    pwln=log(pw)
    pw5=wv^0.5
    presc=press/1013.25
    
    #coefficients for Ebn
    ai11=0.000261031
    ai12=-0.001810295
    ai21=-0.000334819
    ai22=-0.001533816
    ai31=-0.000118343
    ai32=0.000100560
    bi11=0.019908235
    bi12=0.195865498
    bi21=0.019649806
    bi22=0.030513908
    bi31=0.003246244
    bi32=0.010001152
    ci11=0.035765704
    ci12=0.737572097
    ci21=0.017074308
    ci22=0.021712639
    ci31=0.004844209
    ci32=0.005917312
    di11=0.072875947
    di12=1.003594313
    di21=0.009537179
    di22=-0.000136749
    di31=0.001883135
    di32=0.000378418
    
    #coefficients for tauG
    ag11=0.000193733
    ag12=-0.002501003
    ag21=-0.000156219
    ag22=-0.000277670
    ag31=-0.000059040
    ag32=0.000006851
    bg11=-0.001159980
    bg12=0.058524158
    bg21=0.002548229
    bg22=0.004829749
    bg31=0.000788522
    bg32=0.000522673
    cg11=-0.001449466
    cg12=-0.833692155
    cg21=-0.012697614
    cg22=-0.030033430
    cg31=-0.004681306
    cg32=-0.011682439
    dg11=-0.141809777
    dg12=-0.085671381
    dg21=-0.015070590
    dg22=-0.030811277
    dg31=-0.007315151
    dg32=-0.008232404    
    ca1=0.393927007
    ca2=-0.014924316
    ca3=-0.092174236
    ca4=-0.001048172
    ca5=-0.009163093
    ca6=0.006108964
    
    #coefficients for tauB
    ab11=0.000771802
    ab12=-0.007478508
    ab21=-0.000180242
    ab22=-0.000368091
    ab31=-0.000049938
    ab32=0.000049572
    bb11=-0.009565174
    bb12=0.148774990
    bb21=0.002336539
    bb22=0.007087530
    bb31=0.001102680
    bb32=0.000284776
    cb11=0.024325289
    cb12=-1.461925559
    cb21=-0.009964879
    cb22=-0.043313859
    cb31=-0.007992277
    cb32=-0.015281853
    db11=-0.191282122
    db12=-0.099096810
    db21=-0.022005054
    db22=-0.027680914
    db31=-0.006793639
    db32=-0.008818264
    cbb1=0.482261237
    cbb2=-0.016678672
    cbb3=0.914171831
    
    ai1 = ai11 * presc + ai12
    ai2 = ai21 * presc + ai22
    ai3 = ai31 * presc + ai32
    bi1 = bi11 * presc + bi12
    bi2 = bi21 * presc + bi22
    bi3 = bi31 * presc + bi32
    ci1 = ci11 * presc + ci12
    ci2 = ci21 * presc + ci22
    ci3 = ci31 * presc + ci32
    di1 = di11 * presc + di12
    di2 = di21 * presc + di22
    di3 = di31 * presc + di32
       
    ab1 = ab11 * presc + ab12
    ab2 = ab21 * presc + ab22
    ab3 = ab31 * presc + ab32
    bb1 = bb11 * presc + bb12
    bb2 = bb21 * presc + bb22
    bb3 = bb31 * presc + bb32
    cb1 = cb11 * presc + cb12
    cb2 = cb21 * presc + cb22
    cb3 = cb31 * presc + cb32
    db1 = db11 * presc + db12
    db2 = db21 * presc + db22
    db3 = db31 * presc + db32
       
    ag1 = ag11 * presc + ag12
    ag2 = ag21 * presc + ag22
    ag3 = ag31 * presc + ag32
    bg1 = bg11 * presc + bg12
    bg2 = bg21 * presc + bg22
    bg3 = bg31 * presc + bg32
    cg1 = cg11 * presc + cg12
    cg2 = cg21 * presc + cg22
    cg3 = cg31 * presc + cg32
    dg1 = dg11 * presc + dg12
    dg2 = dg21 * presc + dg22
    dg3 = dg31 * presc + dg32
         
    ai=ai1+pw5*ai2+pwln*ai3
    bi=bi1+pw5*bi2+pwln*bi3
    ci=ci1+pw5*ci2+pwln*ci3
    di=di1+pw5*di2+pwln*di3     
      
    ab=ab1+pw5*ab2+pwln*ab3
    bb=bb1+pw5*bb2+pwln*bb3
    cb=cb1+pw5*cb2+pwln*cb3
    db=db1+pw5*db2+pwln*db3

    ag=ag1+pw5*ag2+pwln*ag3
    bg=bg1+pw5*bg2+pwln*bg3
    cg=cg1+pw5*cg2+pwln*cg3
    dg=dg1+pw5*dg2+pwln*dg3
    
    ET=ai*aod3 + bi*aod2 + ci*aod + di
    taub=(ab * aod3 + bb * aod2 + cb * aod + db)
    taug=(ag * aod3 + bg * aod2 + cg * aod + dg)    
    taub=-taub  
    taug=-taug 
          
    b = cbb1*(pw^cbb2)*(cbb3^aod)
    g = ca1 + ca2*pwln + ca3*aodln + ca4*pwln*pwln + ca5*aodln*aodln + ca6*pwln*aodln
    ETcor=ET*En0
    
    #direct normal radiation
    Ebnsolis18=ETcor*exp(-taub/((cos(sza))^b))
    
    #global horizontal irradiance
    Eghsolis18=ETcor*exp(-taug/((cos(sza))^g))*(cos(sza))
	
	#diffuse horizontal radiation    
    Edhsolis18=Eghsolis18-Ebnsolis18*(cos(sza))
    
    ###Quality control
    lower=0
    Ebnsolis18[Ebnsolis18<lower]=0
    Edhsolis18[Edhsolis18<lower]=0
    Eghsolis18[Eghsolis18<lower]=0
	return(list(Ebnsolis18, Edhsolis18, Eghsolis18))
	}


