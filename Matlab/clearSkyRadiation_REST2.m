%% Test example
% number_of_sites = 500000;
% zenith_angle = repmat(deg2rad(-90:90)',[1,number_of_sites]);
% pressure = repmat(linspace(960,1012,size(zenith_angle,1))',[1,number_of_sites]);
% water_vapour = repmat(linspace(7,2,size(zenith_angle,1))',[1,number_of_sites]);
% ozone = repmat(linspace(0.02,0.06,size(zenith_angle,1))',[1,number_of_sites]);
% nitrogen_dioxide = repmat(linspace(0.0002,0.0003,size(zenith_angle,1))',[1,number_of_sites]);
% AOD550 = repmat(linspace(0.2,0.4,size(zenith_angle,1))',[1,number_of_sites]);
% Angstrom_exponent = repmat(linspace(1.1,1.3,size(zenith_angle,1))',[1,number_of_sites]);
% surface_albedo = repmat(linspace(0.3,0.35,size(zenith_angle,1))',[1,number_of_sites]);
% datevecs = datevec(linspace(datenum('01012018','ddmmyyyy'),datenum('02012018','ddmmyyyy'),size(zenith_angle,1)));
% Esc=1366.1;
% ndd = datenum(datevecs) - datenum([datevecs(:,1),ones(length(datevecs(:,1)),2),zeros(length(datevecs(:,1)),3)]);
% beta=(2.*pi.*ndd)./365;
% Eext=Esc*(1.00011+0.034221*cos(beta)+0.00128*sin(beta)+0.000719*cos(2*beta)+0.000077*sin(2*beta));
% Eext = repmat(Eext,[1,number_of_sites]);
% [ghi, dni, dhi] = clearSkyRadiation_REST2(zenith_angle, Eext, pressure, water_vapour, ozone, nitrogen_dioxide, AOD550, Angstrom_exponent, surface_albedo);
% clearvars -except ghi dni dhi

function [ghi, dni, dhi] = clearSkyRadiation_REST2(zenith_angle, Eext, pressure, water_vapour, ozone, nitrogen_dioxide, AOD550, Angstrom_exponent, surface_albedo)
%% NEW CLEAR-SKY MODEL
% This model is called the REST2v5 clear sky model, written and designed by
% Christian A. Gueymard over a series of publications, though primarily in
% his 2008 paper in the Journal of Solar Energy (volume 82, issue 3, pages 
% 272-285) titled "REST2: High-performance solar radiation model for 
% cloudless-sky irradiance,  illuminance, and photosynthetically active 
% radiation - Validation with a benchmark dataset"

% This model has a wide variety of inputs. The temporal operation will
% alter where these come from. They are staged as required inputs instead
% of any decision making occuring within this function. That said, the
% input types and locations will be discussed here.
%
%% Input requirements for the REST2 
% Expected input type is double (though could work for single).
%   zenith_angle        [radians]
%   surface_albedo      [fraction]
%   pressure            [mb]                (local barometric) 
%   Angstrom_exponent   [dimensionless]     (also known as alpha)
%   AOD_550             [dimensionless]     (aerosol optical depth at 550 nm)
%   ozone               [atm.cm]            (total columular amount) 
%   nitrogen_dioxide    [atm.cm]            (total columular amount) 
%   water_vapour        [atm.cm]            (total columular amount) 
%   Eext                [Wm-2]              (Extraterrestrial irradiance)
%
%% Outputs from the REST2
%   ghi (double) global horizontal irradiance  , W m^-2
%   dni (double) direct normal irradiance      , W m^-2
%   dhi (double) diffuse horizontal irradiance , W m^-2
%
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Sources of data %%%%%%%%%%%%%%%%%%%%%%%%%
%% MERRA-2
% Historically, all of these with exception of nitrogen_dioxide are
% available back to 1980 from the MERRA-2 database with 1-hour resolution.
% --------------------------- MERRA-2 Links ------------------------------
% https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/docs/
% https://disc.gsfc.nasa.gov/datasets?keywords=%22MERRA-2%22&page=1&subject=Aerosols&source=Models%2FAnalyses%20MERRA-2&temporalResolution=1%20hour
% https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf

% ---------------------------  MERRA-2 Data ------------------------------
% Variable              Collection  Code        Conversion   Units 
% 
% total_aerosol_ext     M2T1NXAER  TOTSCATTAU   -            - 
% AOD_550               M2T1NXAER  TOTEXTTAU    -            -
% Angstrom_exponent     M2T1NXAER  TOTANGSTR    -            -
% surface_albedo        M2T1NXRAD  ALBEDO       -            frac. 
% ozone                 M2T1NXSLV  TO3          *0.001       atm-cm
% water_vapour          M2T1NXSLV  TQV          *0.1         atm-cm 
% Pressure              M2T1NXSLV  PS           *0.01        hPa or mb
%
% The collecton code is MERRA-2 (M2), Temporal resolution 1-h (T1). Native
% horizontal (NX), and the grouping attributed (AER - aerosol mixing ratio,
% RAD - radiation, SLV - single level.)
% MERRA-2 data files are provided in netCDF-4 format. 
%
%% Reanalysis height correction
% There is need for altitude corrections to the MERRA-2 gridded data due to
% the natural decrease of AOD_559, total_aerosol_ext and water_vapour with 
% increasing altitude.
% This is particularly significant when the reported reanalysis gridded 
% cell elevation differs from the ground station elevation height. 
% The proposed correction was defined by Gueymard (2009) and is based on a 
% scale-height approximation:

% k(h)=k(h_0 ) exp[(h_0-h)/H_a ]

% where k(h) is the variable k at surface height h, k(h_0) is the variable at
% MERRA-2 cell height h_0, and H_a is the scale height found to be suitable
% at a value of 2100m. 
% Gueymard, C.A., Thevenard, D. 2009. Monthly average clear-sky broadband 
% irradiance database for worldwide solar heat gain and building cooling 
% load calculations. Solar Energy. 83, 1998-2018.
% The cell height h_0 can be taken from the constants of MERRA-2:
% the v5.12.4 2-dimensional constants dataset by extracting the surface
% geopotential height (PHIS variable) and dividing it by the average 
% standard gravity (9.80665 ms$^{-2}) resulting in a lookup of h_0. 
% a columular lat-lon-h0 has been provided in a comma delimited .txt file
% called MERRA2-cell-height.txt.
% This cell height correction MUST be applied to the MERRA-2 AOD550, Tau
% and water_vapour BEFORE passing it into this function.
% Recommended usage would be to take the lat lons of the locations and
% perform a knnsearch with the lat lons of the cell height. It is this
% delta that allows for a correction.

%% OMI Nitrogen
% The only variable not extractable from MERRA-2 is the nitrogen dioxide
% amount. This is, however, available from NASA's OMI product that is found
% aboard the Aura polar orbiting satellite from the OMNO2d product, and is 
% converted from OMI units [mol/cm^2] to a columnar amount in 
% atmosphere centimetres [atm-cm], through a division of 2.69*10^16. 
% The OMI data source has a daily resolution and is quality controlled.
% As the OMI instrument is on-board a polar orbiting satellite, there are
% naturally absent recordings in-between flight paths. 
% REST2 only works for nitrogen_dioxide between 0 and 0.03 atm-cm.
% It is possible that a default value of 0.0002	atm-cm is used in absence
% of data.

%% The REST2v5 Model

%% limitations of input
Angstrom_exponent(Angstrom_exponent>2.5)=2.5;
Angstrom_exponent(Angstrom_exponent<0)=0;
pressure(pressure>1100)=1100;
pressure(pressure<300)=300;
water_vapour(water_vapour>10)=10;
water_vapour(water_vapour<0)=0;
ozone(ozone>0.6)=0.6;
ozone(ozone<0)=0;
nitrogen_dioxide(nitrogen_dioxide>0.03)=0.03;
nitrogen_dioxide(nitrogen_dioxide<0)=0;
surface_albedo(surface_albedo>1)=1;
surface_albedo(surface_albedo<0)=0;

%air mass for aerosols extinction
ama=abs((cos(zenith_angle)+0.16851.*(zenith_angle.*180./pi).^0.18198./(95.318-zenith_angle.*180./pi).^1.9542).^-1);
%air mass for water vapor absorption
amw=abs((cos(zenith_angle)+0.10648.*(zenith_angle.*180./pi).^0.11423./(93.781-zenith_angle.*180./pi).^1.9203).^-1);
%air mass for nitrogen dioxide absorption
% amn=abs((cos(zenith_angle)+1.1212.*(zenith_angle.*180./pi).^1.6132./(111.55-zenith_angle.*180./pi).^3.2629).^-1);
%air mass for ozone absorption
amo=abs((cos(zenith_angle)+1.0651.*(zenith_angle.*180./pi).^0.6379./(101.8-zenith_angle.*180./pi).^2.2694).^-1);
%air mass for Rayleigh scattering and uniformly mixed gases absorption
amR=abs((cos(zenith_angle)+0.48353.*((zenith_angle.*180./pi).^0.095846)./(96.741-zenith_angle.*180./pi).^1.754).^-1);
amRe=abs((pressure./1013.25).*(cos(zenith_angle)+0.48353.*((zenith_angle.*180./pi).^0.095846)./(96.741-zenith_angle.*180./pi).^1.754).^-1);

% Angstrom turbidity
ang_beta=AOD550./(0.55.^-Angstrom_exponent);
ang_beta(ang_beta>1.1) = 1.1;
ang_beta(ang_beta<0)=0;

%% Band 1 %%
%transmittance for Rayleigh scattering
TR1=(1 + 1.8169.*amRe - 0.033454.*amRe.^2 )./(1 + 2.063.*amRe + 0.31978.*amRe.^2);
%transmittance for uniformly mixed gases absorption
Tg1=(1 + 0.95885.*amRe + 0.012871.*amRe.^2) ./(1 + 0.96321.*amRe + 0.015455.*amRe.^2);
%transmittance for Ozone absorption
uo=ozone;
f1=uo.*(10.979 - 8.5421.*uo)./(1 + 2.0115.*uo + 40.189.*uo.^2);
f2=uo.*(-0.027589 - 0.005138.*uo)./(1 - 2.4857.*uo + 13.942.*uo.^2);
f3=uo.*(10.995 - 5.5001.*uo)./(1 + 1.6784.*uo + 42.406.*uo.^2);
To1=(1 + f1.*amo + f2.*amo.^2)./(1 + f3.*amo);
%transmittance for Nitrogen dioxide absorption
un=nitrogen_dioxide;
g1=(0.17499 + 41.654.*un - 2146.4.*un.^2)./(1 + 22295.*un.^2);
g2=un.*(-1.2134 + 59.324.*un)./(1 + 8847.8.*un.^2);
g3=(0.17499 + 61.658.*un + 9196.4.*un.^2)./(1 + 74109.*un.^2);
Tn1=min(1,((1 + g1.*amw +  g2.*amw.^2)./(1 + g3.*amw)));
Tn1166=min(1,((1 + g1.*1.66 +  g2.*1.66.^2)./(1 + g3.*1.66)))  ; %at air mass=1.66
%transmittance for Water Vapor absorption
h1=water_vapour.*(0.065445 + 0.00029901.*water_vapour)./(1 + 1.2728.*water_vapour);
h2=water_vapour.*(0.065687 + 0.0013218.*water_vapour)./(1 + 1.2008.*water_vapour);
Tw1=(1 + h1.*amw)./(1 + h2.*amw);
Tw1166=(1 + h1.*1.66)./(1 + h2.*1.66) ;    %at air mass=1.66

%coefficients of angstrom_alpha
AB1=ang_beta;
alph1=Angstrom_exponent;
d0=0.57664 - 0.024743.*alph1;
d1=(0.093942 - 0.2269.*alph1 + 0.12848.*alph1.^2)./(1 + 0.6418.*alph1);
d2=(-0.093819 + 0.36668.*alph1 - 0.12775.*alph1.^2)./(1 - 0.11651.*alph1);
d3=alph1.*(0.15232 - 0.087214.*alph1 + 0.012664.*alph1.^2)./(1 - 0.90454.*alph1 + 0.26167.*alph1.^2);
ua1=log(1 + ama.*AB1);
lam1=(d0 + d1.*ua1 + d2.*ua1.^2)./(1 + d3.*ua1.^2);

%Aeroso transmittance
ta1=abs(AB1.*lam1.^-alph1);
TA1=exp(-ama.*ta1);
%Aeroso scattering transmittance
TAS1=exp(-ama.*0.92.*ta1);    %w1=0.92 recommended

%forward scattering fractions for Rayleigh extinction
BR1=0.5.*(0.89013 - 0.0049558.*amR + 0.000045721.*amR.^2);

%Aerosol scattering correction factor
g0=(3.715 + 0.368.*ama + 0.036294.*ama.^2)./(1 + 0.0009391.*ama.^2);
g1=(-0.164 - 0.72567.*ama + 0.20701.*ama.^2)./(1 + 0.0019012.*ama.^2);
g2=(-0.052288 + 0.31902.*ama + 0.17871.*ama.^2)./(1 + 0.0069592.*ama.^2);
F1=(g0 + g1.*ta1)./(1 + g2.*ta1);

%sky albedo
rs1=(0.13363 + 0.00077358.*alph1 + AB1.*(0.37567 + 0.22946.*alph1)./(1-0.10832.*alph1))./(1 + AB1.*(0.84057 + 0.68683.*alph1)./(1 - 0.08158.*alph1));
%ground albedo
rg=surface_albedo;

%% Band 2 %%
%transmittance for Rayleigh scattering
TR2=(1 - 0.010394.*amRe)./(1-0.00011042.*amRe.^2);
%transmittance for uniformly mixed gases absorption
Tg2=(1 + 0.27284.*amRe - 0.00063699.*amRe.^2)./(1 + 0.30306.*amRe);
%transmittance for Ozone absorption
To2=1; % Ozone (none)
%transmittance for Nitrogen dioxide absorption
Tn2=1 ;% Nitrogen (none)
Tn2166=1;     %at air mass=1.66
%transmittance for water vapor  absorption
c1=water_vapour.*(19.566 - 1.6506.*water_vapour + 1.0672.*water_vapour.^2)./(1 + 5.4248.*water_vapour + 1.6005.*water_vapour.^2);
c2=water_vapour.*(0.50158 - 0.14732.*water_vapour + 0.047584.*water_vapour.^2)./(1 + 1.1811.*water_vapour + 1.0699.*water_vapour.^2);
c3=water_vapour.*(21.286 - 0.39232.*water_vapour + 1.2692.*water_vapour.^2)./(1 + 4.8318.*water_vapour + 1.412.*water_vapour.^2);
c4=water_vapour.*(0.70992 - 0.23155.*water_vapour + 0.096514.*water_vapour.^2)./(1 + 0.44907.*water_vapour + 0.75425.*water_vapour.^2);
Tw2=(1 + c1.*amw + c2.*amw.^2)./(1 + c3.*amw + c4.*amw.^2);
Tw2166=(1 + c1.*1.66 + c2.*1.66.^2)./(1 + c3.*1.66 + c4.*1.66.^2);

%coefficients of angstrom_alpha
AB2=ang_beta;
alph2=Angstrom_exponent;
e0=(1.183 - 0.022989.*alph2 + 0.020829.*alph2.^2)./(1 + 0.11133.*alph2);
e1=(-0.50003 - 0.18329.*alph2 + 0.23835.*alph2.^2)./(1 + 1.6756.*alph2);
e2=(-0.50001 + 1.1414.*alph2 + 0.0083589.*alph2.^2)./(1 + 11.168.*alph2);
e3=(-0.70003 - 0.73587.*alph2 + 0.51509.*alph2.^2)./(1 + 4.7665.*alph2);
ua2 = log(1 + ama.*AB2);
lam2 = (e0 + e1.*ua2 + e2.*ua2.^2)./(1 + e3.*ua2);

%Aeroso transmittance
ta2=abs(AB2.*lam2.^(-alph2));
TA2=exp(-ama.*ta2);
TAS2=exp(-ama.*0.84.*ta2) ;  %w2=0.84 recommended

%forward scattering fractions for Rayleigh extinction
BR2=0.5;  % multi scatter negibile in Band 2
%the aerosol forward scatterance factor
Ba=1 - exp(-0.6931-1.8326.*cos(zenith_angle));

%Aerosol scattering correction
h0=(3.4352 + 0.65267.*ama + 0.00034328.*ama.^2)./(1 + 0.034388.*ama.^1.5);
h1=(1.231 - 1.63853.*ama + 0.20667.*ama.^2)./(1 + 0.1451.*ama.^1.5);
h2=(0.8889 - 0.55063.*ama + 0.50152.*ama.^2)./(1 + 0.14865.*ama.^1.5);
F2=(h0 + h1.*ta2)./(1 + h2.*ta2);

%sky albedo
rs2=(0.010191 + 0.00085547.*alph2 + AB2.*(0.14618 + 0.062758.*alph2)./(1 - 0.19402.*alph2))./(1 + AB2.*(0.58101 + 0.17426.*alph2)./(1 - 0.17586.*alph2));

%irradiance BAND1
E0n1=Eext.*0.46512;
%direct beam irradiance
Ebn1=E0n1.*TR1.*Tg1.*To1.*Tn1.*Tw1.*TA1;
%the incident diffuse irradiance on a perfectly absorbing ground
Edp1=E0n1.*cos(zenith_angle).*To1.*Tg1.*Tn1166.*Tw1166.*(BR1.*(1-TR1).*TA1.^0.25 + Ba.*F1.*TR1.*(1-TAS1.^0.25));
%multiple reflections between the ground and the atmosphere
Edd1=rg.*rs1.*(Ebn1.*cos(zenith_angle) + Edp1)./(1-rg.*rs1);

%irradiance BAND2
E0n2 = Eext.*0.51951;
%direct beam irradiance
Ebn2 = E0n2.*TR2.*Tg2.*To2.*Tn2.*Tw2.*TA2;
%the incident diffuse irradiance on a perfectly absorbing ground
Edp2 = E0n2.*cos(zenith_angle).*To2.*Tg2.*Tn2166.*Tw2166.*(BR2.*(1-TR2).*TA2.^0.25 + Ba.*F2.*TR2.*(1-TAS2.^0.25));
%multiple reflections between the ground and the atmosphere
Edd2 = rg.*rs2.*(Ebn2.*cos(zenith_angle) + Edp2)./(1-rg.*rs2);

%TOTALS BAND1+BAND2
%direct horizontal irradiance
Ebh=(Ebn1 + Ebn2) .* cos(zenith_angle);
dni=Ebn1 + Ebn2;
%correct for zenith angle 
dni(rad2deg(zenith_angle)>90)=0;
Ebh(rad2deg(zenith_angle)>90)=0;

%diffuse horizontal irradiance
dhi=Edp1 + Edd1 + Edp2 + Edd2;
dhi(rad2deg(zenith_angle)>90)=0;

%global horizontal irradiance
ghi = Ebh + dhi;

%Quality Control
lower=0;
ghi(ghi<lower)=NaN;
dni(dni<lower)=NaN;
dhi(dhi<lower)=NaN;

end