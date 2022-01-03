import numpy as np


def clear_sky_mac2(sza, Earth_radius, pressure, wv, ang_beta, ang_alpha, albedo, components):
    """clear_sky_model mac2 1982"""
    # components = 1, output = Edn
    # components = 2, output = [Edn, Edh]
    # components = 3, output = [Egh, Edn, Edh]

    # Extraterrestrial irradiance
    Esc = 1353  # author set 1353
    Eext = Esc * np.power(Earth_radius, -2)

    # Air Mass
    amm = 35 / np.power((1224 * np.power(np.cos(np.deg2rad(sza)), 2) + 1), 0.5)
    amm[amm < 0] = 0

    # Ozone Transmittance
    ozone = 0.35  # Davies and Mckay 1982 set ozone a fixed value of 3.5mm
    XO = amm * (ozone * 10)  # Davies and Mckay 1982 ozone unit is mm, here in the code unit is cm
    AO = ((0.1082 * XO) / (np.power((1 + 13.86 * XO), 0.805))) + ((0.00658 * XO) / (1 + np.power((10.36 * XO), 3))) + (
            0.002118 * XO / (1 + 0.0042 * XO + 3.23e-6 * np.power(XO, 2)))
    TO = 1 - AO

    # Rayleigh Transmittance  (linear interpolation based on  Table 2 in Davies and Mckay 1982 )
    TR = amm
    amms = np.array([0.5, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6, 10, 30]).reshape(16, 1)
    TRs = np.array(
        [.9385, .8973, .8830, .8696, .8572, .8455, .8344, .7872, .7673, .7493, .7328, .7177, .7037, .6907, .6108,
         .4364]).reshape(16, 1)
    TR = np.interp(amm[:, 0], amms[:, 0], TRs[:, 0])
    TR = TR.reshape(np.size(sza, 0), 1)

    # Aerosols Transmittance, borrowed from BH 1981
    tA = ang_beta * (np.power(0.38, -ang_alpha) * 0.2758 + np.power(0.5, -ang_alpha) * 0.35)
    TA = np.exp(-tA * amm)

    # Water Vapor Transmittance
    XW = amm * wv * 10 * np.power((pressure / 1013.25), 0.75)
    AW = 0.29 * XW / (np.power((1 + 14.15 * XW), 0.635) + 0.5925 * XW)
    # Forward Scatter
    szajiao = sza
    szajiaos = np.array([0, 25.8, 36.9, 45.6, 53.1, 60.0, 66.4, 72.5, 78.5, 90]).reshape(10, 1)
    fs = np.array([.92, .91, .89, .86, .83, .78, .71, .67, .60, .60]).reshape(10, 1)
    f = np.interp(szajiao[:, 0], szajiaos[:, 0], fs[:, 0])
    f = f.reshape(np.size(sza, 0), 1)

    lower = 0
    if components == 1:
        # Direct normal irradiance
        EbnMAC2 = Eext * (TO * TR - AW) * TA
        EbnMAC2[EbnMAC2 < lower] = lower  # Quality control

        output = EbnMAC2
    elif components == 2:
        # Direct normal irradiance
        EbnMAC2 = Eext * (TO * TR - AW) * TA
        EbnMAC2[EbnMAC2 < lower] = lower  # Quality control

        # Diffuse horizontal irradiance
        # diffuse components from Rayleigh scatter
        DR = Eext * np.cos(np.deg2rad(sza)) * TO * (1 - TR) / 2
        # diffuse components from scattering by aerosol
        DA = Eext * np.cos(np.deg2rad(sza)) * (TO * TR - AW) * (
                1 - TA) * 0.75 * f  # w0 = 0.75 according to Table5 in Davies and Mckay 1982
        # diffuse horizontal irradiance
        Taaa = np.power(0.95,
                        1.66)  # Taaa is TA determined at amm=1.66, k=0.95 according to Table5 in Davies and Mckay 1982
        poub = 0.0685 + (1 - Taaa) * 0.75 * (
                1 - 0.83)  # f' is f determined at amm=1.66, f' equals  0.83, estimate theta when amm=1.66, theta near 53 degree
        EdhMAC2 = poub * albedo * (EbnMAC2 * np.cos(np.deg2rad(sza)) + DR + DA) / (1 - poub * albedo) + DR + DA
        EdhMAC2[EdhMAC2 < lower] = lower;  # Quality control

        output = [EbnMAC2, EdhMAC2]
    else:
        # Direct normal irradiance
        EbnMAC2 = Eext * (TO * TR - AW) * TA
        EbnMAC2[EbnMAC2 < lower] = lower  # Quality control

        # Diffuse horizontal irradiance
        # diffuse components from Rayleigh scatter
        DR = Eext * np.cos(np.deg2rad(sza)) * TO * (1 - TR) / 2
        # diffuse components from scattering by aerosol
        DA = Eext * np.cos(np.deg2rad(sza)) * (TO * TR - AW) * (
                1 - TA) * 0.75 * f  # w0 = 0.75 according to Table5 in Davies and Mckay 1982
        # diffuse horizontal irradiance
        Taaa = np.power(0.95,
                        1.66)  # Taaa is TA determined at amm=1.66, k=0.95 according to Table5 in Davies and Mckay 1982
        poub = 0.0685 + (1 - Taaa) * 0.75 * (
                1 - 0.83)  # f' is f determined at amm=1.66, f' equals  0.83, estimate theta when amm=1.66, theta near 53 degree
        EdhMAC2 = poub * albedo * (EbnMAC2 * np.cos(np.deg2rad(sza)) + DR + DA) / (1 - poub * albedo) + DR + DA
        EdhMAC2[EdhMAC2 < lower] = lower;  # Quality control

        # Global horizontal irradiance
        EghMAC2 = (EbnMAC2 * np.cos(np.deg2rad(sza)) + DR + DA) / (1 - poub * albedo)
        EghMAC2[EghMAC2 < lower] = lower  # Quality control
        output = [EghMAC2, EbnMAC2, EdhMAC2]

    return output

def model_test():
    size_sza = 181
    sza = np.arange(-90, 91).reshape([size_sza, 1])
    Earth_radius = np.linspace(1, 1.1, size_sza).reshape([size_sza, 1])
    pressure = np.linspace(960, 1012, size_sza).reshape([size_sza, 1])
    wv = np.linspace(7, 2, size_sza).reshape([size_sza, 1])
    ang_beta = np.linspace(1.1, 1.3, size_sza).reshape([size_sza, 1])
    ang_alpha = np.linspace(1, 1.1, size_sza).reshape([size_sza, 1])
    albedo = np.linspace(0.3, 0.35, size_sza).reshape([size_sza, 1])
    components = 3
    output = clear_sky_mac2(sza, Earth_radius, pressure, wv, ang_beta, ang_alpha, albedo, components)
    return output


if __name__ == '__main__':
    print(model_test())
