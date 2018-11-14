%% NEW CLEAR-SKY MODEL
% This model had the best performing GHI from a comprehensive study of 57
% different clear-sky models. It was devised by Davies and Mckay in 1981 at
% the University of MacMaster, Canada. The MMAC model was published in the
% Journal of Solar Energy in a paper titled "Estimating solar irradiance
% and components" 1982, volume 29, issue 1, pages 55-64.
% Further adaptation was made by Gueymard 2003 in his clear-sky validation
% paper where considerations for Aerosol Transition were improved.
%
% This model has a wide variety of inputs. The temporal operation will
% alter where these come from. They are staged as required inputs instead
% of any decision making occuring within this function. That said, the
% input types and locations will be discussed here.
%
%% Input requirements for the MAC model
% Expected input type is double (though could work for single).
%   zenith_angle        [radians]
%   surface_albedo      [fraction]
%   pressure            [mb]                (local barometric)
%   Angstrom_exponent   [dimensionless]     (also known as alpha)
%   AOD_550             [dimensionless]     (aerosol optical depth at 550 nm)
%   water_vapour        [atm.cm]            (total columular amount)
%   Eext                [Wm-2]              (Extraterrestrial irradiance)
% IMPORTANT NOTE ON Eext, THE MAC MODEL USES A SOLAR CONSTANT OF 1353 Wm-2
%
%% Outputs from the REST2
%   ghi (double) global horizontal irradiance  , W m.^-2
%   dni (double) direct normal irradiance      , W m.^-2
%   dhi (double) diffuse horizontal irradiance , W m.^-2
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Sources of data %%%%%%%%%%%%%%%%%%%%%%%%%
%% MERRA-2
% Historically, all of these with exception of nitrogen_dioxide are
% available back to 1980 from the MERRA-2 database with 1-hour resolution.
% --------------------------- MERRA-2 Links ------------------------------
% https:././gmao.gsfc.nasa.gov./reanalysis./MERRA-2./docs./
% https:././disc.gsfc.nasa.gov./datasets?keywords=%22MERRA-2%22&page=1&subject=Aerosols&source=Models%2FAnalyses%20MERRA-2&temporalResolution=1%20hour
% https:././gmao.gsfc.nasa.gov./pubs./docs./Bosilovich785.pdf
%
% ---------------------------  MERRA-2 Data ------------------------------
% Variable              Collection  Code        Conversion   Units
%
% total_aerosol_ext     M2T1NXAER  TOTSCATTAU   -            -
% AOD_550               M2T1NXAER  TOTEXTTAU    -            -
% Angstrom_exponent     M2T1NXAER  TOTANGSTR    -            -
% surface_albedo        M2T1NXRAD  ALBEDO       -            frac.
% water_vapour          M2T1NXSLV  TQV          .*0.1         atm-cm
% Pressure              M2T1NXSLV  PS           .*0.01        hPa or mb
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
%
% k(h)=k(h_0 ) exp[(h_0-h)./H_a ]
%
% where k(h) is the variable k at surface height h, k(h_0) is the variable at
% MERRA-2 cell height h_0, and H_a is the scale height found to be suitable
% at a value of 2100m.
% Gueymard, C.A., Thevenard, D. 2009. Monthly average clear-sky broadband
% irradiance database for worldwide solar heat gain and building cooling
% load calculations. Solar Energy. 83, 1998-2018.
% The cell height h_0 can be taken from the constants of MERRA-2:
% the v5.12.4 2-dimensional constants dataset by extracting the surface
% geopotential height (PHIS variable) and dividing it by the average
% standard gravity (9.80665 ms$.^{-2}) resulting in a lookup of h_0.
% a columular lat-lon-h0 has been provided in a comma delimited .txt file
% called MERRA2-cell-height.txt.
% This cell height correction MUST be applied to the MERRA-2 AOD550, Tau
% and water_vapour BEFORE passing it into this function.
% Recommended usage would be to take the lat lons of the locations and
% perform a knnsearch with the lat lons of the cell height. It is this
% delta that allows for a correction.
%
%% Test example (requirement to have latlon2solzrzenazi)
% lat = -90:10:90;
% leg = cell(size(lat));
% for i = 1:length(lat);
%     leg{i} = num2str(lat(i));
% end
% lon = zeros(size(lat));
% number_of_sites = length(lat);
% datevecs = datevec(linspace(datenum('01012018','ddmmyyyy'),datenum('02012018','ddmmyyyy'),1440));
% Esc=1353;
% ndd = datenum(datevecs) - datenum([datevecs(:,1),ones(length(datevecs(:,1)),2),zeros(length(datevecs(:,1)),3)]);
% beta=(2.*pi.*ndd)./365;
% Eext=Esc.*(1.00011+0.034221.*cos(beta)+0.00128.*sin(beta)+0.000719.*cos(2.*beta)+0.000077.*sin(2.*beta));
% Eext = repmat(Eext,[1,number_of_sites]);
% zenith_angle = deg2rad(latlon2solarzenazi(lat,lon,datevecs));
% % arbitrary input variables.
% pressure = repmat(linspace(960,1012,size(zenith_angle,1))',[1,number_of_sites]);
% water_vapour = repmat(linspace(7,2,size(zenith_angle,1))',[1,number_of_sites]);
% AOD550 = repmat(linspace(0.2,0.4,size(zenith_angle,1))',[1,number_of_sites]);
% Angstrom_exponent = repmat(linspace(1.1,1.3,size(zenith_angle,1))',[1,number_of_sites]);
% surface_albedo = repmat(linspace(0.3,0.35,size(zenith_angle,1))',[1,number_of_sites]);
% [ghi, dni, dhi] = clearSkyRadiation_MAC(zenith_angle, Eext, pressure, water_vapour, AOD550, Angstrom_exponent, surface_albedo);
% % make a figure showing the GHI DNI and DHI 
% figure('Name','Example plot of GHI, DNI and DHI from the MAC model','color','w')
% subplot(1,3,1)
% plot(datetime(datevecs),ghi)
% ylim([0 1200])
% title('GHI')
% ylabel('Irradiance [Wm^-^2]')
% axis square
% subplot(1,3,2)
% plot(datetime(datevecs),dni)
% title('DNI')
% ylim([0 1200])
% axis square
% subplot(1,3,3)
% plot(datetime(datevecs),dhi)
% title('DHI')
% ylim([0 1200])
% legend(leg)
% axis square
% % % expand the axis for better view

function [ghi, dni, dhi] = clearSkyRadiation_MAC(zenith_angle, Eext, pressure, water_vapour, AOD550, Angstrom_exponent, surface_albedo)
% Angstrom turbidity
Angstrom_turbidity = AOD550 ./ (0.55.^-Angstrom_exponent); % not in MAC table but a derived input

%relative optical Air Mass
amm = 35 ./ ( (1224 .* cos(zenith_angle).^2 + 1 ).^0.5 ); % table 2. m_r
amm(amm<0) = NaN;

%Ozone Transmittance. table 2 - transmittance after ozone absporption. 
XO = amm.*3.5; % X1 in table 2.   %Davies and Mckay use a fixed value of 3.5mm, it could potentially improve with dynamic input but then may be overfit
TO = 1 - (( (0.1082.*XO) ./ ((1+13.86.*XO).^0.805) ) - ( (0.00658.*XO) ./ (1+(10.36.*XO).^3) ) - (0.002118.*XO ./ (1+0.0042.*XO+3.23e-6.*XO.^2) ));

%Rayleigh Transmittance . Table 2 - Raylegh scatter R. 
% this is done with lookup tables in the paper, however, we empoy linear
% interpolation of the values within the lookup table as follows.
amms = [0.5, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6, 10, 30];
TRs = [.9385, .8973, .8830, .8696, .8572, .8455, .8344, .7872, .7673, .7493, .7328, .7177, .7037, .6907, .6108, .4364];
TR = interp1(amms,TRs,amm); 

%Aerosols Transmittance from Gueymard 2003 eq (3): Direct solar
%transmittance and irradiance with broadband models. Part I: detailed
%theoretical performance assessment
tA = Angstrom_turbidity .* (0.38.^-Angstrom_exponent .* 0.2758 + 0.5.^-Angstrom_exponent .* 0.35);
TA = exp(-tA .* amm);

% MAC version
% eq 9 in MAC paper: TA = exp(-tA*amm) = k^amm
% values of k and w_o are defined in table 5 as
% 
% STATION  |   k  | 
% goose    | 0.97 | 
% charlotte| 0.95 | etc.

%Water Vapor Transmittance
XW = amm .* water_vapour .* 10 .* (pressure ./ 1013.25) .^ 0.75;
AW = 0.29 .* XW ./ ( ( (1 + 14.15 .* XW) .^ 0.635 ) + 0.5925.*XW);
% table 2 has it using the dewpoint temperature
% UW = exp(2.2572 + 0.05454 .* Td).* (pressure./1013.25).^(0.75).*(To/T).^1/2;
% X2 = amm.* UW
% AW = 0.29 .* X2./(1+14.15.*X2).^0.635 + 0.5925.*X2

%Forward Scatter. In table 2 as ratio of toward to total scatter for
%aerosol
% the forward scatter operates the same as the loop up tables previously.
zenith_lookup = [0, 25.8, 36.9, 45.6, 53.1, 60.0, 66.4, 72.5, 78.5, 90];
fs = [.92, .91, .89, .86, .83, .78, .71, .67, .60, .60];
f = interp1(zenith_lookup, fs, rad2deg(zenith_angle));

%direct beam irradiance. Eq 6
dni = Eext.*(TO.*TR-AW).*TA;

%diffuse components from Rayleigh scatter. Eq 7
DR = Eext.*cos(zenith_angle).*TO.*(1-TR)./2; 
%diffuse components from scattering by aerosol. Eq 8 with 0.75 assumption
w0 = 0.75;
DA = Eext.*cos(zenith_angle).*(TO.*TR-AW).*(1-TA).*w0.*f ;    % 0.75 value for w0

%global horizontal irradiance
TRa = 0.95.^1.66; % Tr(a) in eq 14.
alphab_R=0.0685; % after eq 13.
f_prime = 0.83; %0.83 is the value when f approaches m = 1.6
CA = 1; % cloud amount = 1 at clear skies
alphac_bar = 1; % cloud scattering = 1 in clear skies.
alphab_A = alphab_R+(1-TRa).*w0.*(1-f_prime); % from eq 13. 
alphab = alphab_R.*(1-CA) + alphab_A + alphac_bar.*CA;
ghi = (dni.*cos(zenith_angle) + DR + DA) ./ (1 - alphab.*surface_albedo); % albedo = alpha_s in eq 10

% diffuse horizontal is stated to be the difference between global and
% direct beam irradiance. 
% This is a space for debate as it depends where the surface reflection
% gets considered. 
% options for Diffuse are
dhi = ghi - dni.*cos(zenith_angle); % this will assume reflected irradiance is in the diffuse
%dhi = ( DR + DA) ./ (1 - poub.*surface_albedo); % we lose here the reflected contribution to dni.
%dhi = DR + DA; % this leaves a discrepancy of reflected irradiance.

%Quality Control
lower = 0;
ghi(ghi<lower) = NaN;
dni(dni<lower) = NaN;
dhi(dhi<lower) = NaN;
% The MAC does not produce values of 0 after sunset and can be infact >0
ghi(rad2deg(zenith_angle)>90)=NaN;
dni(rad2deg(zenith_angle)>90)=NaN;
dhi(rad2deg(zenith_angle)>90)=NaN;

end