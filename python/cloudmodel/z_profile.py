import numpy as np

from .constants import (
    Av0,
    Cp,
    Cv,
    G,
    Kapa,
    Lsl0,
    Lvl0,
    P00,
    Rd,
    Rv,
    T0,
    Tmax,
    Tmin,
    Vis0,
    Vtnie0,
    celcius_temperature_aux,
    dx1,
    k,
    nz1,
    rhogra,
)


def latent_heat() -> list[np.ndarray]:
    """
    Calculo de los calores latentes de evaporacion, fusion y sublimacion
    calor latente de sublimacion, Pruppacher and Klett, 2010, eq 4-74
    los coeficientes no estan referenciados
    Inputs:
    Tmin: temperatura minima, en C, int
    Tmax: temperatura maxima, en C, int
    """
    tlvl = Lvl0 * (T0 / k) ** (0.167 + 3.67e-4 * k)

    tlsl = (
        Lsl0
        + 0.485 * celcius_temperature_aux
        - 2.5e-3 * celcius_temperature_aux**2.0
    ) * 4180.0

    tlvs = tlvl + tlsl

    return [tlvl, tlsl, tlvs]


def saturated_vapor_pressure2(
    tlvl: np.array, tlvs: np.array
) -> list[np.ndarray]:
    """
      MAL!!! DA NEGATIVO AUN CON LOS CALCULOS A MANO. MAL PRUPPACHER!!!
      Coeficientes para la tension de vapor de saturacion liquido y solido vs T
      Ceficientes y Expresiones extraidas de Pruppacher and Klett, 2010,
      Apendix A4-8
      Liquido vapor. Validos para -50C < T < 50C. Coeficiente redondeados a
      float. a0 en mB
      Para T< -50C, correcciones usando Clausius-Clayperon, Pruppacher and
      Klett, 2010, eq 4-86
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

    telvs = auxiliar_fun(lv_coefficients)
    tesvs = auxiliar_fun(sv_coefficients)
    for k_temp, i in enumerate(range(10), start=210):
        aux_telvs = tlvl[10] / Rv * (1.0 / 220.0 - 1.0 / k_temp)
        telvs[i] = telvs[10] * np.exp(aux_telvs)
        aux_tesvs = tlvs[10] / Rv * (1.0 / 220.0 - 1.0 / k_temp)
        tesvs[i] = tesvs[10] * np.exp(aux_tesvs)
    for i in range(10):
        telvs[i] = 0.0
        tesvs[i] = 0.0
    return [telvs, tesvs]


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

    tvis = 4.9e-8 * celcius_temperature_aux + Vis0

    tvis[k < T0] -= 1.2e-10 * celcius_temperature_aux[k < T0] ** 2.0
    return tvis


def crystal_efficiencies() -> list[np.ndarray]:
    # Se usa en el rango T_C < T0

    eautcn = 10.0 ** (celcius_temperature_aux * 0.035 - 0.7)
    eacrcn = np.exp(celcius_temperature_aux * 0.09)

    return [eautcn, eacrcn]


def velocities() -> list[np.ndarray]:
    u_z_initial = np.zeros(nz1 + 7)
    v_z_initial = np.zeros(nz1 + 7)

    for i in range(3, nz1 + 4):
        z_aux = i * dx1
        if z_aux <= 500.0:
            u_z_initial[i] = 0.0
            v_z_initial[i] = 0.0
        elif z_aux <= 2000.0:
            z_reference = z_aux - 500.0
            aux = 4.0 * (z_reference / 1500.0) ** 2
            u_z_initial[i] = aux
        elif z_aux <= 9000.0:
            z_reference = z_aux - 2000.0
            base_horizontal_velocity = z_reference / 7000
            u_z_initial[i] = 4.0 - 10.0 * base_horizontal_velocity**2
            v_z_initial[i] = 3.0 * np.sqrt(base_horizontal_velocity)
        else:
            z_reference = z_aux - 9000.0
            u_z_initial[i] = 4.0 * (z_reference / 9000.0) ** 2.0 - 6.0
            v_z_initial[i] = 3.0 - 5.0 * np.sqrt(z_reference / 9000.0)

    u_z_initial = u_z_initial * 0.7
    v_z_initial = u_z_initial * 0.7
    return [u_z_initial, v_z_initial]


def tt_f(z_aux) -> np.array:
    """
    Interpola valores de las temperaturas entre capas, en forma cuadratica
    La derivada primera de la temperatura es lineal y a traves de ella se
    calcula T
    Condicionalmente inestable (lineal) hasta los 2000 y estable entre los
    5500m y los 9000m (lineal)
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
def pp():  # sourcery skip: inline-immediately-returned-variable
    """
    Calcula la presion del aire seco no perturbada
    Se basa en la ecuacion de equilibrio hidrostatico integrando G/(Rd * T(z))
    La integracion se hace sobre un dominio mas fino (4 veces) y para mas
    altura (que deberia ser superflua, ie con nx4 = nz1*4 deberia andar)
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
    for i in range(1, nx4 + 2):
        zetaa = (2 * i - 2) * dx4
        zetam = (2 * i - 1) * dx4
        zetad = (2 * i) * dx4

        ya = 1 / tt_f(zetaa)
        ym = 1 / tt_f(zetam)
        yd = 1 / tt_f(zetad)
        integ[i] = integ[i - 1] + ya + 4 * ym + yd

    presi0 = np.zeros(nz1 + 7)
    for i in range(len(presi0) - 3):
        presi0[i + 3] = P00 * np.exp(-G / Rd * (integ[2 * i] * dx4 / 3))
    for i in range(3):
        presi0[i] = P00
    return presi0


def temperature():  # sourcery skip: inline-immediately-returned-variable
    temperature_z_initial = np.zeros(nz1 + 7)
    for i in range(len(temperature_z_initial) - 3):
        temperature_z_initial[i + 3] = tt_f(i * dx1)
    for i in range(3):
        temperature_z_initial[i] = P00
    return temperature_z_initial


def air_density(presi0, temperature_z_initial):
    # sourcery skip: inline-immediately-returned-variable
    air_density_z_initial = np.zeros(nz1 + 7)
    for i in range(len(air_density_z_initial) - 1):
        air_density_z_initial[i] = presi0[i] / Rd / temperature_z_initial[i]
    for i in range(3):
        air_density_z_initial[i] = presi0[2] / Rd / temperature_z_initial[2]
    return air_density_z_initial


def aerosol():  # sourcery skip: inline-immediately-returned-variable
    aerosol_z_initial = np.zeros(nz1 + 7)
    for i in range(len(aerosol_z_initial) - 3):
        aerosol_z_initial[i + 3] = 10000.0 * np.exp(-(i * dx1) / 2500.0)
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


def vapor(temperature_z_initial, telvs):
    """
    Calcula la densidad de vapor de agua para todos los niveles z
    Telvs: Tension de vapor liquido vapor saturada, kg/m**3 , float array
    (en temperaturas)
    temp: perfil de temperaturas, en K, float array
    rel: humedad relativa, [0, 1], float array, en z
    Rv: Constante de los gases para el vapor de agua
    """
    vapor_z_initial = np.zeros(nz1 + 7)
    for i in range(3, nz1 + 4):
        z_aux = (i - 3) * dx1
        temperature_aux = temperature_z_initial[i]
        n = int(temperature_aux)
        aux = temperature_aux - n
        sat_press_lv_aux = (
            telvs[n - 210] * (1 - aux) + telvs[n + 1 - 210] * aux
        )
        relative_humidity_aux = humidity(z_aux)
        vapor_z_initial[i] = (
            relative_humidity_aux * sat_press_lv_aux / Rv / temperature_aux
        )
    return vapor_z_initial


def air_density_recalc(air_density_z_initial, vapor_z_initial):
    air_density_z_initial_recalc = np.zeros(nz1 + 7)
    for i in range(3):
        air_density_z_initial_recalc[i] = air_density_z_initial[i]
    for i in range(3, nz1 + 1):
        air_density_z_initial_recalc[i] = (
            air_density_z_initial[i] + vapor_z_initial[i]
        )
    for i in range(nz1 + 1, nz1 + 7):
        air_density_z_initial_recalc[i] = air_density_z_initial[i]
    return air_density_z_initial_recalc


def rain_terminal_velocity(presi0):
    """
    Velocidad terminal para gota de lluvia, cte que depende de P
    se define para niveles intermedios
    """
    av = np.zeros(2 * nz1 + 8)
    for i in range(2, nz1 + 3):
        av[2 * i - 1] = (
            Av0
            * ((P00 / presi0[i + 1]) ** 0.286 + (P00 / presi0[i + 2]) ** 0.286)
            / 2.0
        )
        av[2 * i] = Av0 * (P00 / presi0[i + 2]) ** 0.286
    av = np.insert(av, 0, 0)
    return av


def snow_terminal_velocity(presi0):
    """
    Velocidad terminal para la nieve, cte que depende de P
    se define para niveles intermedios
    """
    vtnie = np.zeros(2 * nz1 + 8)
    for i in range(2, nz1 + 3):
        vtnie[2 * i - 1] = (
            Vtnie0
            * ((P00 / presi0[i + 1]) ** 0.3 + (P00 / presi0[i + 2]) ** 0.3)
            / 2.0
        )
        vtnie[2 * i] = Vtnie0 * (P00 / presi0[i + 2]) ** 0.3
    vtnie = np.insert(vtnie, 0, 0)
    return vtnie


def hail_terminal_velocity(tvis, temperature_z_initial, air_density_z_initial):
    """
    Velocidad terminal para el granizo, depende de la viscosidad y de la
    densidad se define para niveles intermedios
    """

    # viscosidad se define por temperaturas y Den por alturas!!
    vtgra0 = np.zeros(2 * nz1 + 8)
    for i in range(2, nz1 + 3):
        aux = 2.754 * rhogra**0.605
        vtgra0[2 * i] = (
            aux
            / tvis[int(temperature_z_initial[i + 2]) - 210] ** 0.21
            / air_density_z_initial[i + 2] ** 0.395
        )
    for i in range(3, nz1 + 3):
        vtgra0[2 * i - 1] = (vtgra0[2 * i - 2] + vtgra0[2 * i]) / 2.0
    vtgra0 = np.insert(vtgra0, 0, 0)
    return vtgra0


def pp2(air_density_z_initial, presi0_aux, temperature_z_initial):
    theta_z_initial = np.zeros(nz1 + 7)
    pres00 = np.zeros(nz1 + 7)
    for i in range(2, nz1 + 6):
        theta_z_initial[i] = (
            temperature_z_initial[i] * (P00 / presi0_aux[i]) ** Kapa
        )
        pres00[i] = temperature_z_initial[i] / theta_z_initial[i]

    # --------------------------------------------------------------!!!!!!!!!!!!!
    den00 = np.zeros(3 * nz1 + 7)
    for i in range(3, nz1 + 3):
        den00[2 * i] = air_density_z_initial[i]
        den00[2 * i + 1] = (
            air_density_z_initial[i] + air_density_z_initial[i + 1]
        ) / 2.0
    den00[2 * (nz1 + 3)] = air_density_z_initial[(nz1 + 3)]
    den00[2 * (nz1 + 3) + 1] = (
        2.0 * air_density_z_initial[(nz1 + 3)] - den00[2 * (nz1 + 3) - 1]
    )

    integ = np.zeros(3 * nz1 + 7)
    for i in range(4, nz1 + 4):
        ya = den00[2 * i - 1]
        ym = den00[2 * i]
        yd = den00[2 * i + 1]
        integ[i] = integ[i - 1] + ya + 4 * ym + yd

    presi0 = np.zeros(nz1 + 7)
    for i in range(4, nz1 + 4):
        presi0[i] = P00 - G * integ[i] * dx1 / 6.0
    presi0[2] = P00
    presi0[3] = P00
    # --------------------------------------------------------------¡¡¡¡¡¡¡¡¡¡¡¡¡¡

    cc2 = np.zeros(nz1 + 7)
    for i in range(3, nz1 + 4):
        theta_z_initial[i] = (
            temperature_z_initial[i] * (P00 / presi0[i]) ** Kapa
        )
        pres00[i] = temperature_z_initial[i] / theta_z_initial[i]
        cc2[i] = Cp * Rd * theta_z_initial[i] * pres00[i] / Cv
    return [presi0, theta_z_initial, pres00, cc2]
