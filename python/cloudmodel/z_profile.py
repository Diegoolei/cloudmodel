from .constants import *
import numpy as np
from matplotlib import pyplot as plt
from numpy import asfortranarray


def latent_heat(Tmin, Tmax):
    """
    Calculo de los calores latentes de evaporacion, fusion y sublimacion
    calor latente de sublimacion, Pruppacher and Klett, 2010, eq 4-74
    los coeficientes no estan referenciados
    Inputs:
    Tmin: temperatura minima, en C, int
    Tmax: temperatura maxima, en C, int
    """
    k = np.arange(210, 314, dtype=np.float64)
    print(k)
    celcius_temperature_aux = k - T0
    Tlvl = Lvl0 * (T0 / k) ** (0.167 + 3.67e-4 * k)

    Tlsl = (
        Lsl0
        + 0.485 * celcius_temperature_aux
        - 2.5e-3 * celcius_temperature_aux**2.0
        ) * 4180.0

    Tlvs = Tlvl + Tlsl

    return [asfortranarray(Tlvl), asfortranarray(Tlsl), asfortranarray(Tlvs)]


def saturated_vapor_pressure(Tmin, Tmax):
    """
      Coeficientes y Expresiones extraidas de Pruppacher and Klett, 2010, Apendix A4-2 y A4-3

    Inputs:
      Tmin: temperatura minima, en C, int
      Tmax: temperatura maxima, en C, int
    """
    # Coeficientes: el primero es para agua liquida y el segundo para hielo
    es = np.array([610.70, 610.64])  # presion de vapor para 0 grado
    A = np.array([7.15, 21.88])  # factor para la temperatura en C
    B = np.array([38.25, 7.65])  # Temperatura de minima en K

    TC2K = T0
    T_C = np.arange(Tmin, Tmax)
    tension_saturada = np.zeros([len(T_C), len(es)])

    for i in range(len(es)):
        tension_saturada[:, i] = es[i] * np.exp(
            A[i] * T_C / (T_C + TC2K - B[i])
        )

    return tension_saturada[:, 0], tension_saturada[:, 1]


# @%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%
# ORIGINAL A DESCARTAR
def saturated_vapor_pressure2(Tmin, Tmax, Tlvl, Tlvs):
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
    a = np.array(
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
    b = np.array(
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

    TC2K = T0

    Telvs = np.zeros(Tmax - Tmin)
    Tesvs = np.zeros(Tmax - Tmin)

    # tensiones de vapor de saturacion para -50C < T_C < 50C
    for k in range(-50 - Tmin, Tmax - Tmin):
        T_C = k + Tmin
        aux = 0
        for i in range(5, -1, -1):
            aux = a[i] + T_C * aux
        Telvs[k] = (
            aux * 100.0
        )  # x100 es para pasar a Pascales SM. Cambiar nombre de Telvs a elvs_T

        aux = 0
        for i in range(5, -1, -1):
            aux = b[i] + T_C * aux
        Tesvs[k] = (
            aux * 100.0
        )  # x100 es para pasar a Pascales SM. Cambiar nombre de Telvs a esvs_T

    # tensiones de vapor de saturacion para T_C < -50C
    if Tmin < -50:
        ind = np.arange(-50 - Tmin)
        # ind = np.arange(2)
        aux = (
            Tlvl[-50 - Tmin]
            / Rv
            * (1.0 / (-50 + TC2K) - 1.0 / (-ind - 50 + TC2K))
        )
        Telvs[ind] = Telvs[-50 - Tmin] * np.exp(aux)
        aux = (
            Tlvs[-50 - Tmin]
            / Rv
            * (1.0 / (-50 + TC2K) - 1.0 / (-ind - 50 + TC2K))
        )
        Tesvs[ind] = Tesvs[-50 - Tmin] * np.exp(aux)

    return Telvs, Tesvs


# @%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%@%


def viscosity(Tmin, Tmax):
    T_C = np.arange(Tmin, Tmax)
    TC2K = T0
    Tvis = 4.9e-8 * T_C + Vis0

    # T_C < 0
    Tvis[T_C < TC2K] = Tvis[T_C < TC2K] - 1.2e-10 * T_C[T_C < TC2K] ** 2.0
    return Tvis


def crystal_efficiencies(Tmin, Tmax):
    # Se usa en el rango T_C < T0
    T_C = np.arange(Tmin, Tmax)
    TC2K = T0

    # cambio por las expresiones de Straka
    Eautcn = 10.0 ** (0.035 * (T_C - 0.7))
    Eacrcn = np.exp(0.09 * T_C)

    return Eautcn, Eacrcn


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
    return np.interp(ind2, ind, aux)


def velocidades(zeta, zeta_p, U_p, V_p):
    """
    Interpola valores de las velocidades horizontales, cuadraticamente, entre capas
    zeta: posiciones en metros, float array
    zeta_p: limite de capas metros, float array
    U_p: velocidad del aire en x, en m/s, float array
    V_p: velocidad del aire en y, en m/s, float array
    """
    UU = np.zeros_like(zeta)
    VV = np.zeros_like(zeta)
    for k in range(len(zeta)):
        for j in range(1, len(zeta_p)):
            if zeta_p[j - 1] < zeta[k] <= zeta_p[j]:
                auxU = (
                    (zeta[k] - zeta_p[j - 1]) / (zeta_p[j] - zeta_p[j - 1])
                ) ** 2
                auxV = (
                    (zeta[k] - zeta_p[j - 1]) / (zeta_p[j] - zeta_p[j - 1])
                ) ** 0.5
                UU[k] = (U_p[j] - U_p[j - 1]) * auxU + U_p[j - 1]
                VV[k] = (V_p[j] - V_p[j - 1]) * auxV + V_p[j - 1]
    return UU, VV


def TTT(zeta, zeta_p, TT0, dT_p):  # !perfil de temperaturas no perturbado
    """
    Interpola valores de las temperaturas entre capas, en forma cuadratica
    La derivada primera de la temperatura es lineal y a traves de ella se calcula T
    Condicionalmente inestable (lineal) hasta los 2000 y estable entre los 5500m y los 9000m (lineal)
    zeta: posiciones en metros, float array
    zeta_p: limite de capas metros, float array
    TT0: temperatura del aire en el piso, en K, float
    dT_p TT0: derivada primera de la temperatura del aire, en K/m, float
    """

    dTT = np.interp(zeta, zeta_p, dT_p)
    TT_f = np.zeros_like(zeta)
    TT_f[0] = TT0
    dx = zeta[1] - zeta[0]
    # dx = zeta_p[1] - zeta_p[0]

    # for k in range(1, len(zeta)):
    #   for j in range(len(zeta_p)-1):
    #      if zeta[k]  <=  zeta_p[j]:
    #        TT_f[k] = TT_f[k-1] + dTT[k] * dx
    for k in range(1, len(zeta)):
        TT_f[k] = TT_f[k - 1] + dTT[k] * dx

    return TT_f


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


# Presion para aire seco
def PP(G, Rd, dx, nz1, Pres0, zeta_p, T_0, dT_p):
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
    dx4 = dx / 4.0  # subdivide el espacio entre capas

    zetaa = np.arange(0, nx4 * dx4 - 2 * dx4, dx4)
    zetam = np.arange(dx4, nx4 * dx4 - dx4, dx4)
    zetad = np.arange(2 * dx4, nx4 * dx4, dx4)

    ya = 1 / TTT(zetaa, zeta_p, T_0, dT_p)
    ym = 1 / TTT(zetam, zeta_p, T_0, dT_p)
    yd = 1 / TTT(zetad, zeta_p, T_0, dT_p)

    integ = np.cumsum(ya + 4 * ym + yd) * dx4 / 3  # revisar indices

    return Pres0 * np.exp(-G / Rd * 2 * integ[: 2 * nz1 + 1 : 2])


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


def main():
    T_heat = latent_heat(Tmin, Tmax)
    Tlvl = T_heat[0]
    Tlsl = T_heat[1]
    Tlvs = T_heat[2]
    Telvs, Tesvs = saturated_vapor_pressure2(Tmin, Tmax, Tlvl, Tlvs)
    Tvis = viscosity(Tmin, Tmax)
    Eautcn, Eacrcn = crystal_efficiencies(Tmin, Tmax)
    UU, VV = velocidades(zeta, zeta_p, U_p, V_p)
    TT_f = TTT(zeta, zeta_p, T_0, dT_p)
    Temp0 = np.copy(TT_f)
    Pres = PP(G, Rd, dx, nz1, P00, zeta_p, T_0, dT_p)
    Den0 = Pres / Rd / Temp0
    Tita0 = Temp0 * (P00 / Pres) ** Kapa
    Pres00 = Temp0 / Tita0
    rel1 = humedad(zeta, zeta_p, rel1_p)

    # Calculo del vapor de agua
    Qvap0 = vapor(Telvs, Temp0, Tmin, rel1, Rv)

    # Recalculo de las cantidades base considerando el vapor de agua
    Den0 = Den0 + Qvap0
    Pres2 = PP2(G, dx, Den0, P00)

    Tita0 = Temp0 * (P00 / Pres) ** Kapa
    Pres00 = Temp0 / Tita0

    # revisar asignaciones y los indices
    Av = rain_terminal_velocity(Pres)
    Vtnie = snow_terminal_velocity(Pres)
    Vtgra0 = hail_terminal_velocity(Tvis, Den0)
