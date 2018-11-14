% datevecs = datevec(datenum(2015,01,01):1/1440:datenum(2016,01,01));
% lat = 0;
% lon = 0;
% NO2 = ones(size(datevecs,1),1).*0.0002;
% wv = ones(size(datevecs,1),1).*2;
% ozone = ones(size(datevecs,1),1).*0.3;
% ang_beta = ones(size(datevecs,1),1).*0.9;
% [EghCalinoiu,EbnCalinoiu,EdhCalinoiu] = CalinouClearSky(NO2,wv,ozone,ang_beta,datevecs,lat,lon);
 
function [EghCalinoiu,EbnCalinoiu,EdhCalinoiu] = CalinouClearSky(NO2,wv,ozone,ang_beta,datevecs,lat,lon)

sza = deg2rad(latlon2solarzenazi(lat,lon,datevecs));
sza(rad2deg(sza)>90)=NaN;
Year = datevecs(:,1);
Dayth = datevec2doy(datevecs);
%extraterrestrial irradiance
Esc=1366.1;
if mod(Year,4) == 0
    totaldayofyear = 366;
else 
    totaldayofyear = 365;
end
B=(Dayth-1).*2.*pi./totaldayofyear;
Eext=Esc.*(1.00011+0.034221.*cos(B)+0.00128.*sin(B)+0.000719.*cos(2.*B)+0.000077.*sin(2.*B));
sinalpha = cos(sza);

%Beam normal irradiance
%ozone
mo3=1./( sinalpha+268.45.*((sza./pi.*180).^0.5).*(115.42-sza./pi.*180).^(-3.2922) );
xo3=mo3.*ozone;
To3= (1+8.5951.*xo3+0.2179.*xo3.^2) ./ (1+8.75308.*xo3+0.45.*xo3.^2-0.0004.*xo3.^3);

%nitrogen dioxide
mno2=1./( sinalpha+602.3.*((sza./pi.*180).^0.5).*(117.96-sza./pi.*180).^(-3.4536) );
xno2=mno2.*NO2;
Tno2=0.742+0.258.*exp(-xno2./0.09);

%water vapour
mw=1./( sinalpha+0.031141.*((sza./pi.*180).^0.1).*(92.471-sza./pi.*180).^(-1.3814) );
xw=mw.*wv;
Tw=(1+1.4.*xw.^0.5+0.053.*xw).*(1-0.013.*xw)./(1+1.626414.*xw.^0.5+0.10816267.*xw);

%trace gas absorption
mg=1./( sinalpha+0.45665.*((sza./pi.*180).^0.07).*(96.4836-sza./pi.*180).^(-1.6970) );
xg=mg;
Tg=(1+0.19558.*xg.^0.5)./(1+0.215582.*xg.^0.5+0.0005.*xg);

%Rayleigh scattering
mR=1./( sinalpha+0.45665.*((sza./pi.*180).^0.07).*(96.4836-sza./pi.*180).^(-1.6970) );
xR=mR;
TR=(1+0.1564.*xR+0.0001.*xR.^2)./(1+0.26038.*xR+0.00697.*xR.^2);

%aerosol extinction
ma=1./( sinalpha+0.031141.*((sza./pi.*180).^0.1).*(92.471-sza./pi.*180).^(-1.3814) );
xa=ma.*ang_beta;
Ta=(1-0.046.*xa)./(1+1.73849.*xa+0.79081.*xa.^2);

taub=To3.*Tno2.*Tw.*Tg.*TR.*Ta;
%direct beam irradiance
EbnCalinoiu=Eext.*taub;

%Diffuse irradiance
%aerosol scattering
g=0.5;
f1g05=1.459.*log(1-g)+0.1595.*(log(1-g)).^2+0.4129.*(log(1-g)).^3;
f2g05=0.0783.*log(1-g)-0.3824.*(log(1-g)).^2+0.5874.*(log(1-g)).^3;
gamma_a05=0.92.*(1-0.5.*exp(f1g05.*sinalpha+f2g05.*sinalpha.^2));
g=0.75;
f1g075=1.459.*log(1-g)+0.1595.*(log(1-g)).^2+0.4129.*(log(1-g)).^3;
f2g075=0.0783.*log(1-g)-0.3824.*(log(1-g)).^2+0.5874.*(log(1-g)).^3;
gamma_a075=0.92.*(1-0.5.*exp(f1g075.*sinalpha+f2g075.*sinalpha.^2));
gamma_a=0.5.*(gamma_a05+gamma_a075);

omega_hat=0.95;
Tas=(0.96865-1.57909.*log(omega_hat)-0.07119.*xa+2.78513.*(log(omega_hat)).^2+0.00157.*xa.^2+0.07939.*xa.*log(omega_hat))./(1-1.22067.*log(omega_hat)+1.83599.*xa+3.20558.*(log(omega_hat)).^2+0.55453.*xa.^2-1.73402.*xa.*log(omega_hat));
Taa=Ta./Tas.*(0.14578+0.69109./omega_hat);

xo3_1=-1+2./0.3.*(xo3-0.15);
mo3_1=-1+2./17.*(mo3-0.01);
cij=[...
    0.907, -0.031, 1.026E-3, -1.106E-4, 1.703E-5, -4.738E-6, 4.409E-6, 1.145E-6, -1.49E-6, 1.748E-8;...
    -0.044, -0.014, 8.96E-4, -8.85E-5, 5.172E-6, -3.739E-6, -7.45E-7, -3.646E-6, 1.707E-6, 7.464E-8;...
    0.023, 7.475E-3, -3.653e-4, 3.38E-5, -1.041e-6, -4.724e-6, -8.628e-7, 2.768e-7, 2.305e-6, 0;...
    -0.012, -4.014E-3, 1.514E-4, -1.806E-5, 2.288E-6, -2.995E-6, -3.292E-6, -1.753E-6, 0, 0;...
    6.22E-3, 2.046E-3, -5.799E-5, 8.207E-6, -1.759E-6, -1.807E-6, 8.908E-7, 0, 0, 0;...
    3.158E-3, -9.618E-4, 2.104E-5, -4.049E-6, 3.605E-6, -1.828E-6, 0, 0, 0, 0;...
    1.576E-3, 4.2E-4, -1.083E-5, 7.229E-6, -3.751E-9, 0, 0, 0, 0, 0;...
    -8.093E-4, -1.768E-4, 6.057E-6, -1.74E-6, 0, 0, 0, 0, 0, 0;...
    4.428E-4, 6.802E-5, -4.685E-6, 0, 0, 0, 0, 0, 0, 0;...
    -2.212E-4, -3.167E-5, 0, 0, 0, 0, 0, 0, 0, 0];

Td_o3=0;
for i = 0:9
    for j = 0:(9-i)
        Xjxo3=cos(j.*acos(xo3_1));
        Ximo3=cos(i.*acos(mo3_1));
        Td_o3=Td_o3+cij((i+1),(j+1)).*Xjxo3.*Ximo3;
    end
end

tauda=gamma_a.*Td_o3.*Tno2.*Tw.*Tg.*TR.*Taa.*(1-Tas);

%Rayleigh
taudR=0.5.*Td_o3.*Tno2.*Tw.*Tg.*Taa.*(1-TR);
%diffuse horizontal irradiance
EdhCalinoiu=Eext.*sinalpha.*(taudR+tauda);

%Global irradiance
EghCalinoiu=EbnCalinoiu.*sinalpha+EdhCalinoiu;

%%%Quality control
lower=0;
EbnCalinoiu(EbnCalinoiu <lower)=0;
EdhCalinoiu(EdhCalinoiu<lower)=0;
EghCalinoiu(EghCalinoiu<lower)=0;

end


function day_of_year=datevec2doy(datevector)
%% safety tests
[rows,cols]=size(datevector);

if cols~=6
    error('datevectors should be produced with datevec() function. This delivers a strict 6 column format. Ensure that the input vectors are as expected.')
end

%% the function
begining_of_year=datenum([datevector(:,1),ones(rows,1),ones(rows,1),zeros(rows,1),zeros(rows,1),zeros(rows,1)]); %derive the beginning of the queried datevec's year
time_now=datenum([datevector(:,1:3),zeros(rows,1),zeros(rows,1),zeros(rows,1)]); %determine the datenum of the queried datevec
day_of_year=time_now-begining_of_year+1;   %Calculate the day of year from that datevec and store

end

function [zen, azi, ha_d] = latlon2solarzenazi(lat, lon, time)
% LATLON2SOLARZENAZI computes zenith & azimuth of sun
% zen = latlon2solarzenazi(lat, lon, time) computes only the zenith, is faster
% INPUTS:
%   lat (array): latitude in degrees (column vector)
%   lon (array): longitude in degrees (column vector)
%   time (datevec): UTC time
% OUTPUTS:
%   zen (array): solar zenith in degrees
%   azi (array): solar azimuth in degrees
%   ha_d (array): hour angle in degrees

if size(lat)~=size(lon)
    error('lat and lon must be same size')
end
if size(lat,1)>size(lat,2)
    warning('lat and lon must be a column vector')
    lat = lat';
    lon = lon';    
end

% broadcast lat and lon to same shape.
% e.g. they may be lat (height x 1) and lon (1 x width)
% `grid = lat .* lon` -> (height x width)
% but `grid(lat > 0)` will only index the first column of the grid!
lat = repmat(lat,[size(time,1),1]);
lon = repmat(lon,[size(time,1),1]);


% UTC timestep: 
% decimal time (e..g 1.5 = 01:30.)
t    = (time(:,4).*60 + time(:,5)); 
dctm = (t./60);
% day of year (this day)
day = datenum( time(:,1:3) ) - datenum( [time(:,1) ones(size(time,1),2)] ) +1;  

% equation of time - minutes
Bd = (360/365.242).*(day-1);                             
Br = deg2rad(Bd);
% equation of time for this day of year, in minutes
eott = (0.258.*cos(Br)-7.416.*sin(Br)-3.648.*cos(2.*Br)-9.228.*sin(2.*Br)); 
% local solar noon in UTC time
usn = 12 - lon./15 -1.*(eott./60) ;                                            

% hour angle
ha_d = (dctm-usn).*15;  
ha_d(ha_d>=180)  = -1 .* (360-ha_d(ha_d>=180));
ha_d(ha_d<=-180) = 360 + ha_d(ha_d<=-180);
ha_r = deg2rad(ha_d);
                                   
% declination
phid_i = (2*pi/365) * (day+(dctm/24)-1);
decl_r = 0.006918...
    - 0.399912*cos(phid_i) + 0.070257*sin(phid_i)...
    - 0.006758*cos(2.*phid_i) + 0.000907*sin(2.*phid_i)...
    - 0.002697*cos(3.*phid_i) + 0.001480*sin(3.*phid_i);

% zenith
tm1 = sin(deg2rad(lat)).*sin(decl_r)+cos(deg2rad(lat)).*cos(decl_r).*cos(ha_r);
zen_r = acos(tm1);
zen = rad2deg(zen_r);

% return only the zenith
if nargout == 1
    return
end

% azimuth (true north=0, positive anti-clockwise)
azi     =   rad2deg(real(asin(sin(ha_r).*cos(decl_r)./sin(zen_r)))); % Duffie Beckman equation 1.6.6b
ha_ew   =   real(acos(tan(decl_r)./tan(deg2rad(lat))));
c1      =   (abs(ha_r) < ha_ew);                          c1(c1==0) = -1;
c1(lat>=0) = c1(lat>=0).*-1;      % THIS SEEMS TO FIX THE BAD NHEM AZI VALS 
c2      =   (deg2rad(lat).*(deg2rad(lat)-decl_r)) >= 0;   c2(c2==0) = -1;
c3      =   (ha_r >= 0);                                  c3(c3==0) = -1;
azi     =   c1.*c2.*azi  +   c3.*((1-c1.*c2)/2).*180;

end


