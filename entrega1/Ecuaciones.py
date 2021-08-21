import numpy
import random

def Drag(v_relativa,magnitud_v_relativa,coeficiente_de_arrastre,var): #todos los datos que se le pasan son datos cuando se esta en t-1
    resultado = -0.75 * (1/var) * coeficiente_de_arrastre * v_relativa * magnitud_v_relativa
    return resultado


def VRelativaX(v_antigua, taus,pz):
    usando = VFluido(taus,pz)
    return v_antigua - usando

def VFluido(taus,pz):
    if 73 * numpy.sqrt(taus) < 5:
        usando = 2.5 * numpy.log(73 * numpy.sqrt(taus) * abs(pz)) + 5.5
    elif 5 <= 73 * numpy.sqrt(taus) < 70:
        usando = 2.5 * numpy.log(73 * numpy.sqrt(taus) * abs(pz)) + 5.5 - 2.5 * numpy.log( 1 + 0.3 * 73 * numpy.sqrt(taus))
    else:
        usando = 2.5 * numpy.log(30 * abs(pz))

    return usando

def MagnitudVRelativa(u_relativa, v_relativa, w_relativa):
    resultado = numpy.sqrt( (u_relativa**2) + (v_relativa**2) + (w_relativa**2) )
    #print(f"MagnitudVRelativa: {resultado:.2f} (u:{u_relativa:.2f},v:{v_relativa:.2f},w:{w_relativa:.2f})")
    # TODO: Sacar condicion de minimo
    #return min(10, resultado)
    return resultado

def CoeficienteDeArrastre(rep):
    componente1 = 1 + (10**4) * (rep**-0.5)
    division1 = 0.208 / componente1
    componente2 = rep * ( 1 + 0.15 * numpy.sqrt(rep) + 0.017 * rep)
    resultado = 24 / (componente2 - division1)
    return resultado

def Rep(magnitud_v_relativa, taus):
    resultado = magnitud_v_relativa * numpy.sqrt(taus) * 73
    return resultado

def PesoSumergidoX(constante, taus, teta):
    resultado = (1/constante) * numpy.sin(teta) * (1/taus)
    return resultado

def PesoSumergidoZ(constante, taus, teta):
    resultado = - (1/constante) * numpy.cos(teta) * (1/taus)
    return resultado

def MasaVirtual(constante, pz_antiguo, w_relativa):
    resultado = (0.5/ constante) *  w_relativa * (2.5/pz_antiguo)
    return resultado

def VSuperior(u_particula,taus,pz, v_relativa, w_relativa):
    u_relativa = u_particula - VFluido(taus,pz+0.5)
    resultado = (u_relativa**2) + (v_relativa**2) + (w_relativa**2)
    #print(f"VSuperior:{resultado:.2f} (u:{u_relativa:.2f},v:{v_relativa:.2f},w:{w_relativa:.2f})")
    return resultado

def VBotton(u_particula,taus,pz, v_relativa, w_relativa):
    u_relativa = u_particula - VFluido(taus,pz-0.5)
    resultado = (u_relativa**2) + (v_relativa**2) + (w_relativa**2)
    #print(f"VInferior:{resultado:.2f} (u:{u_relativa:.2f},v:{v_relativa:.2f},w:{w_relativa:.2f})")
    return resultado

def FuerzaElevacion(constante,CL, top, bot):
    resultado = 0.75 * (1/constante) * CL * (top - bot)
    #if (top-bot > 20):
    #print(f"FuerzaElevacion: {resultado:.2f} (top:{top:.2f},bot:{bot:.2f},dif:{top-bot:.2f})")
    return resultado

def AnguloRebote(w_prima, u_actual):
    resultado = numpy.arctan(w_prima/u_actual)
    while resultado > numpy.radians(75):
        resultado = numpy.arctan(w_prima/u_actual)
    return resultado

def UPrima(alfa,w_prima):
    error = random.uniform(numpy.radians(0), numpy.radians(10))
    resultado = w_prima / numpy.tan(alfa + error)
    return resultado

def VPrima(u_prima):
    alfa = random.uniform(numpy.radians(-10), numpy.radians(10))
    resultado = u_prima * numpy.tan(alfa)
    return resultado

def Posicion(p_vieja, v_actual, delta_t):
    resultado = p_vieja + v_actual * delta_t
    return resultado
    
def PosicionZ(p_vieja, v_actual, delta_t):
    resultado = p_vieja + v_actual * delta_t
    print(resultado)
    return resultado

def Velocidad(v_vieja, delta_t, fuerzas):
    sumatoria = 0
    for fuerza in fuerzas:
        sumatoria += fuerza
    resultado = v_vieja + delta_t * sumatoria
    return resultado
