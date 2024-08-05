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
    return [Telvs, Tesvs]


def auxiliar_fun(coefficients):
    aux_pre = np.array([0.0] * len(k))
    aux_pre = coefficients[3] + celcius_temperature_aux * (
        coefficients[4]
        + celcius_temperature_aux
        * (coefficients[5] + celcius_temperature_aux * coefficients[6])
    )

    aux = np.array([0.0] * len(k))
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


def TT_f(z_aux: np.array) -> np.array:
    """
    Interpola valores de las temperaturas entre capas, en forma cuadratica
    La derivada primera de la temperatura es lineal y a traves de ella se calcula T
    Condicionalmente inestable (lineal) hasta los 2000 y estable entre los 5500m y los 9000m (lineal)
    zeta: posiciones en metros, float array
    zeta_p: limite de capas metros, float array
    TT0: temperatura del aire en el piso, en K, float
    dT_p TT0: derivada primera de la temperatura del aire, en K/m, float
    """
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

    Presi0 = np.zeros(nz1 + 5)
    for k in range(1, nz1 + 4):
        Presi0[k + 1] = P00 * np.exp(-G / Rd * (integ[2 * k] * dx4 / 3))
    Presi0[0] = P00
    Presi0[1] = P00
    for _ in range(2):
        Presi0 = np.insert(Presi0, 0, P00)
    return Presi0


def temperature():  # sourcery skip: inline-immediately-returned-variable
    temperature_z_initial = np.zeros(nz1 + 5)
    for k in range(-1, nz1 + 4):
        temperature_z_initial[k + 1] = TT_f(k * dx1)
    for _ in range(2):
        temperature_z_initial = np.insert(temperature_z_initial, 0, P00)
    return temperature_z_initial


def air_density(Presi0, temperature_z_initial):
    # sourcery skip: inline-immediately-returned-variable
    air_density_z_initial = np.zeros(nz1 + 7)
    for k in range(len(Presi0) - 4):
        air_density_z_initial[k] = (
            Presi0[k + 4] / Rd / temperature_z_initial[k + 4]
        )
    return air_density_z_initial


def aerosol():  # sourcery skip: inline-immediately-returned-variable
    aerosol_z_initial = np.zeros(nz1 + 7)
    for k in range(len(aerosol_z_initial) - 4):
        aerosol_z_initial[k] = 10000.0 * np.exp(-((k + 1) * dx1) / 2500.0)
    return aerosol_z_initial


# Presion para aire humedo
def PP2(G, dx, Den0, Pres0):
    """
    Calcula la presion del aire humedo no perturbada
    Se basa en la ecuacion de equilibrio hidrostatico integrando G * rho
    La integracion es del tipo Simpson en un dominio mas fino (2 veces)
    G: aceleracion de la gravedad, m/s, float
    dx: espaciamiento de capas en z, m, float
    Den0: densidad del aire humedo en z, kg/m**3, float array
    Pres0: presion a nivel del suelo, Pascales, float
    """

    nz1 = len(Den0)
    ind = np.arange(nz1)
    ind2 = np.arange(nz1 * 2 + 2)  # revisar estos indices (SM 23/7/24)

    Den00 = np.interp(ind2, ind, Den0)

    ya = Den00[:-3:2]
    ym = Den00[1:-2:2]
    yd = Den00[2:-1:2]

    integ = np.cumsum(ya + 4 * ym + yd) * dx / 3.0  # revisar indices

    return Pres0 - G * integ


def rain_terminal_velocity(Presi):  # revisar indices!
    """
    Velocidad terminal para gota de lluvia, cte que depende de P
    se define para niveles intermedios
    """
    nz = len(Presi)
    ind = np.arange(nz)
    ind2 = np.arange(0, nz, 0.5)  # revisar estos indices (SM 23/7/24)
    aux = Av0 * (P00 / Presi) ** 0.286
    return np.interp(ind2, ind, aux)


def snow_terminal_velocity(Presi):
    """
    Velocidad terminal para la nieve, cte que depende de P
    se define para niveles intermedios
    """

    nz = len(Presi)
    ind = np.arange(nz)
    ind2 = np.arange(0, nz, 0.5)  # revisar estos indices (SM 23/7/24)
    aux = Vtnie0 * (P00 / Presi) ** 0.3
    return np.interp(ind2, ind, aux)


def hail_terminal_velocity(Tvis, Den0):
    """
    Velocidad terminal para el granizo, depende de la viscosidad y de la densidad
    se define para niveles intermedios
    """

    ## viscosidad se define por temperaturas y Den por alturas!!

    nz = len(Tvis)
    ind = np.arange(nz)
    ind2 = np.arange(0, nz, 0.5)  # revisar estos indices (SM 23/7/24)
    aux = 2.754 * rhogra**0.605 / Tvis**0.21 / Den0**0.395
    Vtgra0 = np.interp(ind2, ind, aux)
    return Vtgra0


def humedad(zeta, zeta_p, H_p):  # !perfil de humedad relativa no perturbado
    """
    Interpola valores de las humedades relativas entre capas, en forma lineal
    zeta: posiciones en metros, float array
    zeta_p: limite de capas metros, float array
    H_p: humedad relativa, [0, 1], float array
    """
    return np.interp(zeta, zeta_p, H_p)


# ---------------------------------------------------------------
# Lo que sigue no lo puedo chequear pues necesitaria las funciones de Telvs,el valor de Rv, etc
def vapor(Telvs, temp, Tmin, rel, Rv):
    """
    Calcula la densidad de vapor de agua para todos los niveles z
    Telvs: Tension de vapor liquido vapor saturada, kg/m**3 , float array (en temperaturas)
    temp: perfil de temperaturas, en K, float array
    rel: humedad relativa, [0, 1], float array, en z
    Rv: Constante de los gases para el vapor de agua
    """
    Qvap0 = np.zeros_like(temp)
    TC = temp - T0
    for k in range(len(temp)):
        n = int(TC[k] - Tmin)
        aux = TC[k] - Tmin - n
        elv1 = (
            Telvs[n] * (1 - aux) + Telvs[n + 1] * aux
        )  # interpolacion para la tension de vapor saturada liquido vapor
        Qvap0[k] = rel[k] * elv1 / Rv / temp[k]

    return Qvap0


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
    # rel1 = humedad(zeta, zeta_p, rel1_p)

    # Calculo del vapor de agua
    # Qvap0 = vapor(Telvs, Temp0, Tmin, rel1, Rv)

    # Recalculo de las cantidades base considerando el vapor de agua
    # Den0 = Den0 + Qvap0
    # Presi0 = PP2(Den0, Presi0)

    # theta_z_initial = theta(temperature_z_initial, Presi0)
    # Pres00 = Temp0 / Tita0

    # revisar asignaciones y los indices
    # Av = rain_terminal_velocity(Presi0)
    # Vtnie = snow_terminal_velocity(Presi0)
    # Vtgra0 = hail_terminal_velocity(Tvis, Den0)
