from .constants import *
import numpy as np
from matplotlib import pyplot as plt


def latent_heat() -> list[np.ndarray]:
    """
    Calculo de los calores latentes de evaporacion, fusion y sublimacion
    calor latente de sublimacion, Pruppacher and Klett, 2010, eq 4-74
    los coeficientes no estan referenciados
    Inputs:
    Tmin: temperatura minima, en C, int
    Tmax: temperatura maxima, en C, int
    """
    Tlvl = Lvl0 * (T0 / k) ** (0.167 + 3.67e-4 * k)

    Tlsl = (
        Lsl0
        + 0.485 * celcius_temperature_aux
        - 2.5e-3 * celcius_temperature_aux**2.0
    ) * 4180.0

    Tlvs = Tlvl + Tlsl

    return [Tlvl, Tlsl, Tlvs]


def saturated_vapor_pressure2(
    Tlvl: np.array, Tlvs: np.array
) -> list[np.ndarray]:
    # sourcery skip: extract-duplicate-method, inline-immediately-returned-variable
    """
      MAL!!! DA NEGATIVO AUN CON LOS CALCULOS A MANO. MAL PRUPPACHER!!!
      Coeficientes para la tension de vapor de saturacion liquido y solido vs T
      Ceficientes y Expresiones extraidas de Pruppacher and Klett, 2010, Apendix A4-8
      Liquido vapor. Validos para -50C < T < 50C. Coeficiente redondeados a float. a0 en mB
      Para T< -50C, correcciones usando Clausius-Clayperon, Pruppacher and Klett, 2010, eq 4-86
    Inputs:
      Tmin: temperatura minima, en C, int
      Tmax: temperatura maxima, en C, int
    """
    # coeficientes para liquido vapor, a[0] en mB
    lv_coefficients = np.array(
        [
            6.10780,
            4.43652e-1,
            1.42895e-2,
            2.65065e-4,
            3.03124e-6,
            2.03408e-8,
            6.13682e-11,
        ]
    )

    # coeficientes para solido vapor, b[0] en mB
    sv_coefficients = np.array(
        [
            6.10918,
            5.03470e-1,
            1.88601e-2,
            4.17622e-4,
            5.82472e-6,
            4.83880e-8,
            1.83883e-10,
        ]
    )

    Telvs = auxiliar_fun(lv_coefficients)
    Tesvs = auxiliar_fun(sv_coefficients)
    for k_temp, i in enumerate(range(10), start=210):
        aux_Telvs = Tlvl[10] / Rv * (1.0 / 220.0 - 1.0 / k_temp)
        Telvs[i] = Telvs[10] * np.exp(aux_Telvs)
        aux_Tesvs = Tlvs[10] / Rv * (1.0 / 220.0 - 1.0 / k_temp)
        Tesvs[i] = Tesvs[10] * np.exp(aux_Tesvs)
    for i in range(10):
        Telvs[i] = 0.0
        Tesvs[i] = 0.0
    return [Telvs, Tesvs]


def auxiliar_fun(coefficients):
    aux_pre = np.zeros(len(k))
    aux_pre = coefficients[3] + celcius_temperature_aux * (
        coefficients[4]
        + celcius_temperature_aux
        * (coefficients[5] + celcius_temperature_aux * coefficients[6])
    )

    aux = np.zeros(len(k))
    aux = coefficients[0] + celcius_temperature_aux * (
        coefficients[1]
        + celcius_temperature_aux
        * (coefficients[2] + celcius_temperature_aux * aux_pre)
    )
    result = np.zeros(Tmax - Tmin)
    result = aux * 100.0

    return result


def viscosity() -> np.ndarray:

    Tvis = 4.9e-8 * celcius_temperature_aux + Vis0

    Tvis[k < T0] -= 1.2e-10 * celcius_temperature_aux[k < T0] ** 2.0
    return Tvis


def crystal_efficiencies() -> list[np.ndarray]:
    # Se usa en el rango T_C < T0

    Eautcn = 10.0 ** (celcius_temperature_aux * 0.035 - 0.7)
    Eacrcn = np.exp(celcius_temperature_aux * 0.09)

    return [Eautcn, Eacrcn]


def velocities() -> list[np.ndarray]:
    biased_nz1 = nz1 + 4
    u_z_initial = np.zeros(biased_nz1)
    v_z_initial = np.zeros(biased_nz1)

    for k in range(biased_nz1 - nz1 - 1, biased_nz1):
        z_aux = (k - (biased_nz1 - nz1 - 1)) * dx1
        if z_aux <= 500.0:
            u_z_initial[k] = 0.0
            v_z_initial[k] = 0.0
        elif z_aux <= 2000.0:
            z_reference = z_aux - 500.0
            aux = 4.0 * (z_reference / 1500.0) ** 2
            u_z_initial[k] = aux
        elif z_aux <= 9000.0:
            z_reference = z_aux - 2000.0
            base_horizontal_velocity = z_reference / 7000
            u_z_initial[k] = 4.0 - 10.0 * base_horizontal_velocity**2
            v_z_initial[k] = 3.0 * np.sqrt(base_horizontal_velocity)
        else:
            z_reference = z_aux - 9000.0
            u_z_initial[k] = 4.0 * (z_reference / 9000.0) ** 2.0 - 6.0
            v_z_initial[k] = 3.0 - 5.0 * np.sqrt(z_reference / 9000.0)

    u_z_initial = u_z_initial * 0.7
    v_z_initial = u_z_initial * 0.0
    return [u_z_initial, v_z_initial]


def TT_f(z_aux) -> np.array:
    """
    Interpola valores de las temperaturas entre capas, en forma cuadratica
    La derivada primera de la temperatura es lineal y a traves de ella se calcula T
    Condicionalmente inestable (lineal) hasta los 2000 y estable entre los 5500m y los 9000m (lineal)
    zeta: posiciones en metros, float array
    zeta_p: limite de capas metros, float array
    TT0: temperatura del aire en el piso, en K, float
    dT_p TT0: derivada primera de la temperatura del aire, en K/m, float
    """
    np.float32(z_aux)
    a = 298.15
    if z_aux <= 2000:
        return a - 9.0e-3 * z_aux
    elif z_aux <= 5500:
        xx = z_aux - 2000.0
        return a - 18.0 - xx * (9.0e-3 - 2e-3 * xx / 3500.0 / 2.0)
    elif z_aux <= 9000:
        xx = z_aux - 5500.0
        return a - 46.0 - 7e-3 * xx
    elif z_aux <= 11000:
        xx = z_aux - 9000
        return a - 70.5 - 7e-3 * xx + 1.75e-6 * xx**2.0
    elif z_aux <= 12000:
        return a - 77.5
    else:
        xx = z_aux - 12000
        return a - 77.5 + 50.0 * (xx / 9000.0) ** 2.0


# Presion para aire seco
def PP():  # sourcery skip: inline-immediately-returned-variable
    """
    Calcula la presion del aire seco no perturbada
    Se basa en la ecuacion de equilibrio hidrostatico integrando G/(Rd * T(z))
    La integracion se hace sobre un dominio mas fino (4 veces) y para mas altura (que deberia ser superflua, ie con nx4 = nz1*4 deberia andar)
    La integracion es del tipo Simpson
    G: aceleracion de la gravedad, m/s, float
    Rd: constante de los gases ideales para el aire seco, s/(m K), float
    dx: espaciamiento de capas en z, m, float
    nz1: numero de puntos en z
    Pres0: presion a nivel del suelo, Pascales, float
    zeta_p: limite de capas metros, float array
    TT0: temperatura del aire en el piso, en K, float
    dT_p TT0: derivada primera de la temperatura del aire, en K/m, float
    """

    nx4 = 500  # considera 500 niveles en altura para la integracion
    dx4 = dx1 / 4.0  # subdivide el espacio entre capas
    integ = np.zeros(nx4 + 2)
    for k in range(1, nx4 + 2):
        zetaa = (2 * k - 2) * dx4
        zetam = (2 * k - 1) * dx4
        zetad = (2 * k) * dx4

        ya = 1 / TT_f(zetaa)
        ym = 1 / TT_f(zetam)
        yd = 1 / TT_f(zetad)
        integ[k] = integ[k - 1] + ya + 4 * ym + yd

    Presi0 = np.zeros(nz1 + 7)
    for k in range(len(Presi0) - 3):
        Presi0[k + 3] = P00 * np.exp(-G / Rd * (integ[2 * k] * dx4 / 3))
    for i in range(3):
        Presi0[i] = P00
    return Presi0


def temperature():  # sourcery skip: inline-immediately-returned-variable
    temperature_z_initial = np.zeros(nz1 + 7)
    for k in range(len(temperature_z_initial) - 3):
        temperature_z_initial[k + 3] = TT_f(k * dx1)
    for i in range(3):
        temperature_z_initial[i] = P00
    return temperature_z_initial


def air_density(Presi0, temperature_z_initial):
    # sourcery skip: inline-immediately-returned-variable
    air_density_z_initial = np.zeros(nz1 + 7)
    for k in range(len(air_density_z_initial) - 1):
        air_density_z_initial[k] = Presi0[k] / Rd / temperature_z_initial[k]
    for i in range(3):
        air_density_z_initial[i] = Presi0[2] / Rd / temperature_z_initial[2]
    return air_density_z_initial


def aerosol():  # sourcery skip: inline-immediately-returned-variable
    aerosol_z_initial = np.zeros(nz1 + 7)
    for k in range(len(aerosol_z_initial) - 3):
        aerosol_z_initial[k + 3] = 10000.0 * np.exp(-(k * dx1) / 2500.0)
    for i in range(3):
        aerosol_z_initial[i] = -aerosol_z_initial[3]

    return aerosol_z_initial


def humidity(z_aux):
    """
    Interpola valores de las humedades relativas entre capas, en forma lineal
    zeta: posiciones en metros, float array
    zeta_p: limite de capas metros, float array
    H_p: humedad relativa, [0, 1], float array
    """
    if z_aux <= 500:
        relative_humidity_aux = 0.55 + 0.05 * z_aux / 500.0
    elif z_aux <= 1500.0:
        relative_humidity_aux = 0.6
    elif z_aux <= 4000:
        relative_humidity_aux = 0.6 - (z_aux - 1500) / 2500.0 * 0.25
    elif z_aux <= 7000:
        relative_humidity_aux = 0.35 - (z_aux - 4000.0) / 3000.0 * 0.25
    else:
        relative_humidity_aux = 0.1 - (z_aux - 7000) / 3000.0 * 0.02
    return relative_humidity_aux


def vapor(temperature_z_initial, Telvs):
    """
    Calcula la densidad de vapor de agua para todos los niveles z
    Telvs: Tension de vapor liquido vapor saturada, kg/m**3 , float array (en temperaturas)
    temp: perfil de temperaturas, en K, float array
    rel: humedad relativa, [0, 1], float array, en z
    Rv: Constante de los gases para el vapor de agua
    """
    vapor_z_initial = np.zeros(nz1 + 7)
    for k in range(3, nz1 + 4):
        z_aux = (k - 3) * dx1
        temperature_aux = temperature_z_initial[k]
        n = int(temperature_aux)
        aux = temperature_aux - n
        sat_press_lv_aux = Telvs[n - 210] * (1 - aux) + Telvs[n + 1 - 210] * aux
        relative_humidity_aux = humidity(z_aux)
        vapor_z_initial[k] = (
            relative_humidity_aux * sat_press_lv_aux / Rv / temperature_aux
        )
    return vapor_z_initial


def air_density_recalc(air_density_z_initial, vapor_z_initial):
    air_density_z_initial_recalc = np.zeros(nz1 + 7)
    for k in range(len(air_density_z_initial_recalc) - 1):
        air_density_z_initial_recalc[k] = air_density_z_initial(
            k
        ) + vapor_z_initial(k)
    for i in range(3):
        air_density_z_initial_recalc[i] = air_density_z_initial(
            2
        ) + vapor_z_initial(2)
    return air_density_z_initial_recalc


def rain_terminal_velocity(Presi0):  # revisar indices!
    """
    Velocidad terminal para gota de lluvia, cte que depende de P
    se define para niveles intermedios
    """
    Av = np.zeros(2 * nz1 + 8)
    for k in range(2, nz1 + 3):
        Av[2 * k - 1] = (
            Av0
            * ((P00 / Presi0[k+1]) ** 0.286 + (P00 / Presi0[k+2]) ** 0.286)
            / 2.0
        )
        Av[2 * k] = (Av0 * (P00 / Presi0[k+2]) ** 0.286)
    Av = np.insert(Av, 0, 0)
    Av_1 =np.array([
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 1462.21582, 1469.43164,
    1476.78625, 1484.14111, 1491.63843, 1499.13599, 1506.78040, 1514.42480,
    1522.22034, 1530.01599, 1537.96753, 1545.91882, 1554.03052, 1562.14221,
    1570.41846, 1578.69482, 1587.13904, 1595.58337, 1604.19934, 1612.81531,
    1621.60620, 1630.39697, 1639.36646, 1648.33569, 1657.48706, 1666.63855,
    1675.97571, 1685.31287, 1694.83923, 1704.36560, 1714.08484, 1723.80408,
    1733.71985, 1743.63562, 1753.75159, 1763.86768, 1774.18774, 1784.50781,
    1795.03699, 1805.56616, 1816.31042, 1827.05481, 1838.02063, 1848.98621,
    1860.17981, 1871.37354, 1882.80200, 1894.23035, 1905.90039, 1917.57056,
    1929.48987, 1941.40918, 1953.58521, 1965.76099, 1978.20142, 1990.64209,
    2003.35547, 2016.06873, 2029.06348, 2042.05823, 2055.34033, 2068.62231,
    2082.18433, 2095.74634, 2109.57739, 2123.40845, 2137.49536, 2151.58203,
    2165.91016, 2180.23853, 2194.79248, 2209.34619, 2224.10840, 2238.87085,
    2253.83154, 2268.79224, 2283.95312, 2299.11377, 2314.47705, 2329.84033,
    2345.40771, 2360.97485, 2376.74219, 2392.50928, 2408.47119, 2424.43286,
    2440.58325, 2456.73315, 2473.06567, 2489.39795, 2505.90576, 2522.41382,
    0.00000000, 0.00000000, 0.00000000
])
    Av =np.array([
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 1462.21585852, 1469.43171703, 
    1476.78638484, 1484.14105265, 1491.63851749, 1499.13598234, 1506.78038509, 1514.42478784,
    1522.22042987, 1530.0160719, 1537.96742292, 1545.91877394, 1554.03047154, 1562.14216914,
    1570.41846152, 1578.69475391, 1587.13911711, 1595.58348031, 1604.1993574, 1612.8152345,
    1621.60610681, 1630.39697911, 1639.36636402, 1648.33574893, 1657.48719758, 1666.63864622,
    1675.97574092, 1685.31283562, 1694.83918714, 1704.36553865, 1714.08478319, 1723.80402774,
    1733.71982379, 1743.63561985, 1753.75164473, 1763.86766962, 1774.18772196, 1784.5077743,
    1795.03692568, 1805.56607707, 1816.31037164, 1827.05466621, 1838.02038919, 1848.98611216,
    1860.17978851, 1871.37346487, 1882.80187067, 1894.23027648, 1905.90045086, 1917.57062525,
    1929.48988305, 1941.40914084, 1953.58508596, 1965.76103108, 1978.20157066, 1990.64211024,
    2003.35546964, 2016.06882903, 2029.06356766, 2042.0583063, 2055.34022348, 2068.62214067,
    2082.18431054, 2095.74648041, 2109.57746646, 2123.40845251, 2137.49528795, 2151.5821234,
    2165.91036564, 2180.23860789, 2194.79241636, 2209.34622483, 2224.10858414, 2238.87094345,
    2253.83165017, 2268.79235689, 2283.95300612, 2299.11365536, 2314.47691924, 2329.84018313,
    2345.40745439, 2360.97472565, 2376.74203699, 2392.50934832, 2408.47108602, 2424.43282371,
    2440.58302154, 2456.73321937, 2473.06557375, 2489.39792813, 2505.90581505, 2522.41370196,
    0.00000000, 0.00000000, 0.00000000
    ])
    return Av_1

def snow_terminal_velocity(Presi0):
    """
    Velocidad terminal para la nieve, cte que depende de P
    se define para niveles intermedios
    """
    Vtnie = np.zeros(2 * nz1 + 8)
    for k in range(2, nz1 + 3):
        Vtnie[2 * k - 1] = (
            Vtnie0
            * ((P00 / Presi0[k+1]) ** 0.3 + (P00 / Presi0[k+2]) ** 0.3)
            / 2.0
        )
        Vtnie[2 * k] = Vtnie0 * (P00 / Presi0[k+2]) ** 0.3
    Vtnie = np.insert(Vtnie, 0, 0)
    return Vtnie


def hail_terminal_velocity(Tvis, temperature_z_initial, air_density_z_initial):
    """
    Velocidad terminal para el granizo, depende de la viscosidad y de la densidad
    se define para niveles intermedios
    """

    ## viscosidad se define por temperaturas y Den por alturas!!
    Vtgra0 = np.zeros(2 * nz1 + 8)
    for k in range(2, nz1 + 3):
        aux = 2.754 * rhogra**0.605
        Vtgra0[2 * k] = (
            aux
            / Tvis[int(temperature_z_initial[k+2]) - 210] ** 0.21
            / air_density_z_initial[k+2] ** 0.395
        )
    for k in range(3, nz1 + 3):
        Vtgra0[2 * k - 1] = (Vtgra0[2 * k - 2] + Vtgra0[2 * k]) / 2.0
    Vtgra0 = np.insert(Vtgra0, 0, 0)
    return Vtgra0

# Presion para aire humedo
def PP2(air_density_z_initial, Presi0):
    """
    Calcula la presion del aire humedo no perturbada
    Se basa en la ecuacion de equilibrio hidrostatico integrando G * rho
    La integracion es del tipo Simpson en un dominio mas fino (2 veces)
    G: aceleracion de la gravedad, m/s, float
    dx: espaciamiento de capas en z, m, float
    Den0: densidad del aire humedo en z, kg/m**3, float array
    Pres0: presion a nivel del suelo, Pascales, float
    """
    Den00 = np.zeros(2 * nz1 + 2)
    for k in range(nz1 - 1):
        Den00[2 * k] = air_density_z_initial[k]
        Den00[2 * k + 1] = (
            air_density_z_initial(k) + air_density_z_initial(k + 1)
        ) / 2.0
    Den00[2 * nz1] = air_density_z_initial(nz1)
    Den00[2 * nz1 + 1] = 2.0 * air_density_z_initial(nz1) - Den00[2 * nz1 - 1]

    integ = np.zeros(nz1 + 1)
    for k in range(1, nz1):
        ya = Den00[2 * k - 1]
        ym = Den00[2 * k]
        yd = Den00[2 * k + 1]
        integ[k] = integ[k - 1] + ya + 4 * ym + yd

    Pres00 = np.zeros(nz1 + 2)
    for k in range(1, nz1):
        Pres00[k] = Presi0 - G * integ[k] * dx1 / 6.0
    Pres00[0] = Presi0
    Pres00[-1] = Presi0


# ---------------------------------------------------------------
# Lo que sigue no lo puedo chequear pues necesitaria las funciones de Telvs,el valor de Rv, etc


def main():
    T_heat = latent_heat()
    Tlvl = T_heat[0]
    Tlsl = T_heat[1]
    Tlvs = T_heat[2]
    sv_pressure = saturated_vapor_pressure2(Tlvl, Tlvs)
    Telvs = sv_pressure[0]
    Tesvs = sv_pressure[1]
    Tvis = viscosity()
    cry_effic = crystal_efficiencies()
    Eautcn = cry_effic[0]
    Eacrcn = cry_effic[1]
    speed = velocities()
    u_z_initial = speed[0]
    v_z_initial = speed[1]

    Presi0 = PP()
    temperature_z_initial = temperature()
    air_density_z_initial = air_density(Presi0, temperature_z_initial)
    aerosol_z_initial = aerosol()

    vapor_z_initial = vapor(temperature_z_initial, Telvs)

    air_density_z_initial = air_density(air_density_z_initial, vapor_z_initial)

    Av = rain_terminal_velocity(Presi0)
    Vtnie = snow_terminal_velocity(Presi0)
    Vtgra0 = hail_terminal_velocity(
        Tvis, temperature_z_initial, air_density_z_initial
    )

    Presi0 = PP2(air_density_z_initial, Presi0)

    # theta_z_initial = theta(temperature_z_initial, Presi0)
    # Pres00 = Temp0 / Tita0
