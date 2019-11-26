import numpy as np
import math


def clear_sky_reset2(zenith_angle: np.ndarray, Eext: np.ndarray, pressure: np.ndarray, water_vapour: np.ndarray,
                     ozone: np.ndarray, nitrogen_dioxide: np.ndarray, AOD550: np.ndarray, Angstrom_exponent: np.ndarray,
                     surface_albedo: np.ndarray):
    '''Every Variable Need to be np.ndarry. np.matrix will cause fatal error'''

    Angstrom_exponent[Angstrom_exponent > 2.5] = 2.5
    Angstrom_exponent[Angstrom_exponent < 0] = 0
    pressure[pressure > 1100] = 1100
    pressure[pressure < 300] = 300
    water_vapour[water_vapour > 10] = 10
    water_vapour[water_vapour < 0] = 0
    ozone[ozone > 0.6] = 0.6
    ozone[ozone < 0] = 0
    nitrogen_dioxide[nitrogen_dioxide > 0.03] = 0.03
    nitrogen_dioxide[nitrogen_dioxide < 0] = 0
    surface_albedo[surface_albedo > 1] = 1
    surface_albedo[surface_albedo < 0] = 0

    # air mass for aerosols extinction
    complex_temp = np.array(zenith_angle * 180. / np.pi, dtype=np.complex)

    ama = np.abs(np.power(np.cos(zenith_angle) + 0.16851 * np.power(complex_temp, 0.18198) / np.power(
        95.318 - complex_temp, 1.9542), -1))
    # air mass for water vapor absorption
    amw = np.abs(np.power(np.cos(zenith_angle) + 0.10648 * np.power(complex_temp, 0.11423) / np.power(
        93.781 - complex_temp, 1.9203), -1))
    # air mass for nitrogen dioxide absorption
    # amn = np.abs(np.power(np.cos(zenith_angle) + 1.1212 * np.power(zenith_angle * 180. / np.pi, 1.6132) / np.power(
    #   3.2629 - zenith_angle * 180. / np.pi, 1.9203), -1))
    # air mass for ozone absorption
    amo = np.abs(np.power(np.cos(zenith_angle) + 1.0651 * np.power(complex_temp, 0.6379) / np.power(
        101.8 - complex_temp, 2.2694), -1))
    # air mass for Rayleigh scattering and uniformly mixed gases absorption
    amR = np.abs(np.power(np.cos(zenith_angle) + 0.48353 * np.power(complex_temp, 0.095846) / np.power(
        96.741 - complex_temp, 1.754), -1))
    amRe = np.abs((pressure / 1013.25) * np.power(
        np.cos(zenith_angle) + 0.48353 * (np.power(complex_temp, 0.095846)) / np.power(
            96.741 - complex_temp, 1.754), -1))

    # Angstrom turbidity
    ang_beta = AOD550 / np.power(0.55, -1 * Angstrom_exponent)
    ang_beta[ang_beta > 1.1] = 1.1
    ang_beta[ang_beta < 0] = 0

    '''Band 1'''

    # transmittance for Rayleigh scattering
    TR1 = (1 + 1.8169 * amRe - 0.033454 * np.power(amRe, 2)) / (1 + 2.063 * amRe + 0.31978 * np.power(amRe, 2))
    # transmittance for uniformly mixed gases absorption
    Tg1 = (1 + 0.95885 * amRe + 0.012871 * np.power(amRe, 2)) / (1 + 0.96321 * amRe + 0.015455 * np.power(amRe, 2))
    # transmittance for Ozone absorption
    uo = ozone
    f1 = uo * (10.979 - 8.5421 * uo) / (1 + 2.0115 * uo + 40.189 * np.power(uo, 2))
    f2 = uo * (-0.027589 - 0.005138 * uo) / (1 - 2.4857 * uo + 13.942 * np.power(uo, 2))
    f3 = uo * (10.995 - 5.5001 * uo) / (1 + 1.6784 * uo + 42.406 * np.power(uo, 2))
    To1 = (1 + f1 * amo + f2 * np.power(amo, 2)) / (1 + f3 * amo)
    # transmittance for Nitrogen dioxide absorption
    un = nitrogen_dioxide
    g1 = (0.17499 + 41.654 * un - 2146.4 * np.power(un, 2)) / (1 + 22295. * np.power(un, 2))
    g2 = un * (-1.2134 + 59.324 * un) / (1 + 8847.8 * np.power(un, 2))
    g3 = (0.17499 + 61.658 * un + 9196.4 * np.power(un, 2)) / (1 + 74109. * np.power(un, 2))
    Tn1_middle = ((1 + g1 * amw + g2 * np.power(amw, 2)) / (1 + g3 * amw))
    Tn1_middle[Tn1_middle > 1] = 1
    Tn1 = Tn1_middle
    # Tn1 = min(1, ((1 + g1 * amw + g2 * np.power(amw, 2)) / (1 + g3 * amw)))
    Tn1166_middle = (1 + g1 * 1.66 + g2 * np.power(1.66, 2)) / (1 + g3 * 1.66)
    Tn1166_middle[Tn1166_middle > 1] = 1
    Tn1166 = Tn1166_middle
    # Tn1166 = min(1, ((1 + g1 * 1.66 + g2 * np.power(1.66, 2)) / (1 + g3 * 1.66)))  # atairmass = 1.66
    # transmittance for Water Vapor absorption
    h1 = water_vapour * (0.065445 + 0.00029901 * water_vapour) / (1 + 1.2728 * water_vapour)
    h2 = water_vapour * (0.065687 + 0.0013218 * water_vapour) / (1 + 1.2008 * water_vapour)
    Tw1 = (1 + h1 * amw) / (1 + h2 * amw)
    Tw1166 = (1 + h1 * 1.66) / (1 + h2 * 1.66)  # atairmass = 1.66

    # coefficients of angstrom_alpha
    AB1 = ang_beta
    alph1 = Angstrom_exponent
    d0 = 0.57664 - 0.024743 * alph1
    d1 = (0.093942 - 0.2269 * alph1 + 0.12848 * np.power(alph1, 2)) / (1 + 0.6418 * alph1)
    d2 = (-0.093819 + 0.36668 * alph1 - 0.12775 * np.power(alph1, 2)) / (1 - 0.11651 * alph1)
    d3 = alph1 * (0.15232 - 0.087214 * alph1 + 0.012664 * np.power(alph1, 2)) / (
            1 - 0.90454 * alph1 + 0.26167 * np.power(alph1, 2))
    ua1 = np.log(1 + ama * AB1)
    lam1 = (d0 + d1 * ua1 + d2 * np.power(ua1, 2)) / (1 + d3 * np.power(ua1, 2))

    # Aeroso transmittance
    ta1 = np.abs(AB1 * np.power(lam1, -1 * alph1))
    TA1 = np.exp(-ama * ta1)

    # Aeroso scattering transmittance
    TAS1 = np.exp(-ama * 0.92 * ta1)  # w1 = 0.92recommended

    # forward scattering fractions for Rayleigh extinction
    BR1 = 0.5 * (0.89013 - 0.0049558 * amR + 0.000045721 * np.power(amR, 2))

    # Aerosol scattering correction factor
    g0 = (3.715 + 0.368 * ama + 0.036294 * np.power(ama, 2)) / (1 + 0.0009391 * np.power(ama, 2))
    g1 = (-0.164 - 0.72567 * ama + 0.20701 * np.power(ama, 2)) / (1 + 0.0019012 * np.power(ama, 2))
    g2 = (-0.052288 + 0.31902 * ama + 0.17871 * np.power(ama, 2)) / (1 + 0.0069592 * np.power(ama, 2))
    F1 = (g0 + g1 * ta1) / (1 + g2 * ta1)

    # sky albedo
    rs1 = (0.13363 + 0.00077358 * alph1 + AB1 * (0.37567 + 0.22946 * alph1) / (1 - 0.10832 * alph1)) / (
            1 + AB1 * (0.84057 + 0.68683 * alph1) / (1 - 0.08158 * alph1))
    # ground albedo
    rg = surface_albedo

    '''Band 2'''

    # transmittance for Rayleigh scattering
    TR2 = (1 - 0.010394 * amRe) / (1 - 0.00011042 * np.power(amRe, 2))
    # transmittance for uniformly mixed gases absorption
    Tg2 = (1 + 0.27284 * amRe - 0.00063699 * np.power(amRe, 2)) / (1 + 0.30306 * amRe)
    # transmittance for Ozone absorption
    To2 = 1  # Ozone (none)
    # transmittance for Nitrogen dioxide absorption
    Tn2 = 1  # Nitrogen (none)
    Tn2166 = 1  # at air mass=1.66

    # transmittance for water vapor  absorption
    c1 = water_vapour * (19.566 - 1.6506 * water_vapour + 1.0672 * np.power(water_vapour, 2)) / (
            1 + 5.4248 * water_vapour + 1.6005 * np.power(water_vapour, 2))
    c2 = water_vapour * (0.50158 - 0.14732 * water_vapour + 0.047584 * np.power(water_vapour, 2)) / (
            1 + 1.1811 * water_vapour + 1.0699 * np.power(water_vapour, 2))
    c3 = water_vapour * (21.286 - 0.39232 * water_vapour + 1.2692 * np.power(water_vapour, 2)) / (
            1 + 4.8318 * water_vapour + 1.412 * np.power(water_vapour, 2))
    c4 = water_vapour * (0.70992 - 0.23155 * water_vapour + 0.096514 * np.power(water_vapour, 2)) / (
            1 + 0.44907 * water_vapour + 0.75425 * np.power(water_vapour, 2))
    Tw2 = (1 + c1 * amw + c2 * np.power(amw, 2)) / (1 + c3 * amw + c4 * np.power(amw, 2))
    Tw2166 = (1 + c1 * 1.66 + c2 * np.power(1.66, 2)) / (1 + c3 * 1.66 + c4 * np.power(1.66, 2))

    # coefficients of angstrom_alpha
    AB2 = ang_beta
    alph2 = Angstrom_exponent
    e0 = (1.183 - 0.022989 * alph2 + 0.020829 * np.power(alph2, 2)) / (1 + 0.11133 * alph2)
    e1 = (-0.50003 - 0.18329 * alph2 + 0.23835 * np.power(alph2, 2)) / (1 + 1.6756 * alph2)
    e2 = (-0.50001 + 1.1414 * alph2 + 0.0083589 * np.power(alph2, 2)) / (1 + 11.168 * alph2)
    e3 = (-0.70003 - 0.73587 * alph2 + 0.51509 * np.power(alph2, 2)) / (1 + 4.7665 * alph2)
    ua2 = np.log(1 + ama * AB2)
    lam2 = (e0 + e1 * ua2 + e2 * np.power(ua2, 2)) / (1 + e3 * ua2)

    # Aeroso transmittance
    lam2_temp = np.array(lam2, dtype=np.complex)
    ta2 = np.abs(AB2 * np.power(lam2_temp, -1 * alph2))
    TA2 = np.exp(-1 * ama * ta2)
    TAS2 = np.exp(-1 * ama * 0.84 * ta2)  # w2=0.84 recommended

    # forward scattering fractions for Rayleigh extinction
    BR2 = 0.5  # multi scatter negibile in Band 2
    # the aerosol forward scatterance factor
    Ba = 1 - np.exp(-0.6931 - 1.8326 * np.cos(zenith_angle))

    # Aerosol scattering correction
    h0 = (3.4352 + 0.65267 * ama + 0.00034328 * np.power(ama, 2)) / (1 + 0.034388 * np.power(ama, 1.5))
    h1 = (1.231 - 1.63853 * ama + 0.20667 * np.power(ama, 2)) / (1 + 0.1451 * np.power(ama, 1.5))
    h2 = (0.8889 - 0.55063 * ama + 0.50152 * np.power(ama, 2)) / (1 + 0.14865 * np.power(ama, 1.5))
    F2 = (h0 + h1 * ta2) / (1 + h2 * ta2)

    # sky albedo
    rs2 = (0.010191 + 0.00085547 * alph2 + AB2 * (0.14618 + 0.062758 * alph2) / (1 - 0.19402 * alph2)) / (
            1 + AB2 * (0.58101 + 0.17426 * alph2) / (1 - 0.17586 * alph2))

    # irradiance BAND1
    E0n1 = Eext * 0.46512
    # direct beam irradiance
    Ebn1 = E0n1 * TR1 * Tg1 * To1 * Tn1 * Tw1 * TA1
    print(E0n1 * np.cos(zenith_angle) * To1 * Tg1 )
    # the incident diffuse irradiance on a perfectly absorbing ground
    Edp1 = E0n1 * np.cos(zenith_angle) * To1 * Tg1 * Tn1166 * Tw1166 * (
            BR1 * (1 - TR1) * np.power(TA1, 0.25) + Ba * F1 * TR1 * (1 - np.power(TAS1, 0.25)))
    # multiple reflections between the ground and the atmosphere
    Edd1 = rg * rs1 * (Ebn1 * np.cos(zenith_angle) + Edp1) / (1 - rg * rs1)

    # irradiance BAND2
    E0n2 = Eext * 0.51951
    # direct beam irradiance
    Ebn2 = E0n2 * TR2 * Tg2 * To2 * Tn2 * Tw2 * TA2
    # the incident diffuse irradiance on a perfectly absorbing ground
    Edp2 = E0n2 * np.cos(zenith_angle) * To2 * Tg2 * Tn2166 * Tw2166 * (
            BR2 * (1 - TR2) * np.power(TA2, 0.25) + Ba * F2 * TR2 * (1 - np.power(TAS2, 0.25)))
    # multiple reflections between the ground and the atmosphere
    Edd2 = rg * rs2 * (Ebn2 * np.cos(zenith_angle) + Edp2) / (1 - rg * rs2)
    # TOTALS BAND1+BAND2
    # direct horizontal irradiance
    Ebh = (Ebn1 + Ebn2) * np.cos(zenith_angle)
    dni = Ebn1 + Ebn2
    # correct for zenith angle
    dni[np.rad2deg(zenith_angle) > 90] = 0
    Ebh[np.rad2deg(zenith_angle) > 90] = 0
    # diffuse horizontal irradiance
    dhi = Edp1 + Edd1 + Edp2 + Edd2
    dhi[np.rad2deg(zenith_angle) > 90] = 0

    # global horizontal irradiance
    ghi = Ebh + dhi

    # Quality Control
    lower = 0
    ghi[ghi < lower] = np.nan
    dni[dni < lower] = np.nan
    dhi[dhi < lower] = np.nan

    return [ghi, dni, dhi]


def data_eext_builder(number_sites, size_zenith):
    esc = 1366.1
    ndd = np.linspace(0, 1, size_zenith).reshape([size_zenith, 1])  # dayth from 1.1 per year
    beta = (2 * np.pi * ndd) / 365
    Eext = esc * (1.00011 + 0.034221 * np.cos(beta) + 0.00128 * np.sin(beta) + 0.000719 * np.cos(
        2 * beta) + 0.000077 * np.sin(
        2 * beta))
    Eext = np.tile(Eext, number_sites)
    return Eext


def model_test():
    number_sites = 1
    size_zenith = 181
    zenith_angle = np.tile(np.deg2rad(np.arange(-90, 91).reshape([size_zenith, 1])), number_sites)
    pressure = np.tile(np.linspace(960, 1012, size_zenith).reshape([size_zenith, 1]), number_sites)
    water_vapour = np.tile(np.linspace(7, 2, size_zenith).reshape([size_zenith, 1]), number_sites)
    ozone = np.tile(np.linspace(0.02, 0.06, size_zenith).reshape([size_zenith, 1]), number_sites)
    nitrogen_dioxide = np.tile(np.linspace(0.0002, 0.0003, size_zenith).reshape([size_zenith, 1]), number_sites)
    AOD550 = np.tile(np.linspace(0.2, 0.4, size_zenith).reshape([size_zenith, 1]), number_sites)
    Angstrom_exponent = np.tile(np.linspace(1.1, 1.3, size_zenith).reshape([size_zenith, 1]), number_sites)
    surface_albedo = np.tile(np.linspace(0.3, 0.35, size_zenith).reshape([size_zenith, 1]), number_sites)
    Eext = data_eext_builder(number_sites, size_zenith)
    [ghi, dni, dhi] = clear_sky_reset2(zenith_angle, Eext, pressure, water_vapour, ozone, nitrogen_dioxide, AOD550, Angstrom_exponent, surface_albedo)

    print(ghi, dni, dhi)
    return 1
if __name__ == '__main__':
    model_test()
