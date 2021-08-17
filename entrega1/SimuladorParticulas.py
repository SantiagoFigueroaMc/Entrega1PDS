import math
import random
from random import seed
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import time
import threading


# Diccionarios para memorizar
VFLUIDO = {}
REP = {}
MAG_RELATIVAS = {}
PESOSUMERGIDO = {}



class Particula:
    px = 0
    py = 0
    pz = 0
    vx = 0
    vy = 0
    vz = 0

    def CrearParticula(self,p1,p2,p3,v1,v2,v3):
        self.px = float(p1)
        self.py = float(p2)
        self.pz = float(p3)
        self.vx = float(v1)
        self.vy = float(v2)
        self.vz = float(v3)

class Variables:
    delta_t = 0
    tiempo_simulacion = 0
    angulo = 0
    relacion_densidad_agua = 0
    taus = 0
    CL = 0

def Drag(v_relativa,magnitud_v_relativa,coeficiente_de_arrastre,var): #todos los datos que se le pasan son datos cuando se esta en t-1
    #resultado = -0.75 * (1/var) * coeficiente_de_arrastre * v_relativa * magnitud_v_relativa
    return -0.75 * (1/var) * coeficiente_de_arrastre * v_relativa * magnitud_v_relativa


def VRelativaX(v_antigua, taus, pz):
    usando = VFluido(taus,pz)
    return v_antigua - usando


def VFluido(taus,pz):
    v = 0.0
    if (pz in VFLUIDO):
        v = VFLUIDO[pz]
    else:
        if TAUS_TIME_73 < 5:
            v = 2.5 * math.log(TAUS_TIME_73 * abs(pz)) + 5.5
        elif 5 <= TAUS_TIME_73 < 70:
            v = 2.5 * math.log(TAUS_TIME_73 * abs(pz)) + 5.5 - 2.5 * math.log( 1 + 0.3 * TAUS_TIME_73)
        else:
            v = 2.5 * math.log(30 * abs(pz))
        VFLUIDO[pz] = v
    return v

def MagnitudVRelativa(u_relativa, v_relativa, w_relativa):
    #resultado = math.sqrt( (u_relativa**2) + (v_relativa**2) + (w_relativa**2) )
    return math.sqrt( (u_relativa**2) + (v_relativa**2) + (w_relativa**2) )

def CoeficienteDeArrastre(rep):
    componente1 = 1 + (10**4) * (rep**-0.5)
    division1 = 0.208 / componente1
    componente2 = rep * ( 1 + 0.15 * math.sqrt(rep) + 0.017 * rep)
    resultado = 24 / (componente2 - division1)
    return resultado

def Rep(magnitud_v_relativa, taus):
    if (magnitud_v_relativa in MAG_RELATIVAS):
        return MAG_RELATIVAS[magnitud_v_relativa]
    else:
        resultado = magnitud_v_relativa * TAUS_TIME_73
        MAG_RELATIVAS[magnitud_v_relativa] = resultado
        return resultado

def PesoSumergidoX(constante, taus, teta):
    #resultado = (1/constante) * math.sin(teta) * (1/taus)
    return (1/constante) * math.sin(teta) * (1/taus)

def PesoSumergidoZ(constante, taus, teta):
    #resultado = - (1/constante) * math.cos(teta) * (1/taus)
    return - (1/constante) * math.cos(teta) * (1/taus)

def MasaVirtual(constante, pz_antiguo, w_relativa):
    #resultado = (0.5/ constante) *  w_relativa * (2.5/pz_antiguo)
    return (0.5/ constante) *  w_relativa * (2.5/pz_antiguo)

def VSuperior(u_particula, taus, pz, v_relativa, w_relativa):
    u_relativa = u_particula - VFluido(taus,pz+0.5)
    resultado = (u_relativa**2) + (v_relativa**2) + (w_relativa**2)
    return resultado

def VBotton(u_particula, taus, pz, v_relativa, w_relativa):
    u_relativa = u_particula - VFluido(taus,pz-0.5)
    resultado = (u_relativa**2) + (v_relativa**2) + (w_relativa**2)
    return resultado

def FuerzaElevacion(constante, CL, top, bot):
    #resultado = 0.75 * (1/constante) * CL * (top - bot)
    return 0.75 * (1/constante) * CL * (top - bot)

def AnguloRebote(w_prima, u_actual):
    #resultado = math.atan(w_prima/u_actual)
    """
    while resultado > math.radians(75):
        print('help')
        resultado = math.atan(w_prima/u_actual)
    """
    return math.atan(w_prima/u_actual)

def UPrima(alfa,w_prima):
    error = random.uniform(math.radians(0), math.radians(10))
    resultado = w_prima / math.tan(alfa + error)
    return resultado

def VPrima(u_prima):
    alfa = random.uniform(math.radians(-10), math.radians(10))
    resultado = u_prima * math.tan(alfa)
    return resultado

def Posicion(p_vieja, v_actual, delta_t):
    #resultado = p_vieja + v_actual * delta_t
    return p_vieja + v_actual * delta_t

def Velocidad(v_vieja, delta_t, fuerzas):
    sumatoria = sum(fuerzas)
    #for fuerza in fuerzas:
    #     sumatoria += fuerza
    resultado = v_vieja + delta_t * sumatoria
    return resultado

def ObtenerDatos(nombre_archivo):
    contador_linea=0
    particulas = []
    simulacion = Variables()
    archivo = open(nombre_archivo,'r')
    while True:
        linea = archivo.readline()
        linea = linea.split()

        if len(linea) != 0:
            if contador_linea == 0:
                simulacion.tiempo_simulacion = float(linea[0])
                simulacion.delta_t = float(linea[1])

            elif contador_linea == 1:
                simulacion.angulo = float(linea[0])
                simulacion.relacion_densidad_agua = float(linea[1])
                simulacion.taus = float(linea[2])
                simulacion.CL = float(linea[3])

            else:
                particula = Particula()
                particula.CrearParticula(linea[0],linea[1],linea[2],linea[3],linea[4],linea[5])
                particulas.append(particula)

            contador_linea+=1
        else:
            break

    archivo.close()
    return simulacion,particulas

def SimularParticula(p_id, particula_simulada,constante,const, t_0):
    print(f"Simulando particula #{p_id}")

    # Variables necesarias
    tiempo_transcurrido = 0
    total_saltos = 0
    lista_alturas_alcanzadas = []
    h_max = 0.0
    h_promedio = 0.0
    direccion_z = False #True sube, False baja
    direccion_z_antiguo = False
    nueva_particula = particula_simulada

    # Simulacion
    while (tiempo_transcurrido < tiempo_total_simulacion):
        px_antiguo = nueva_particula.px
        py_antiguo = nueva_particula.py
        pz_antiguo = nueva_particula.pz
        vx_antiguo = nueva_particula.vx
        vy_antiguo = nueva_particula.vy
        vz_antiguo = nueva_particula.vz

        # General
        vrX = VRelativaX(vx_antiguo, constante.taus, pz_antiguo)
        vrY = vy_antiguo
        vrZ = vz_antiguo

        mgR = MagnitudVRelativa(vrX, vrY, vrZ)
        rep = Rep(mgR, constante.taus)
        if (rep in REP):
            coeficiente = REP[rep]
        else:
            coeficiente =  CoeficienteDeArrastre(rep)
            REP[rep] = coeficiente

        # Fuerzas X
        fuerza_arrastre_x = Drag(vrX,mgR,coeficiente,const)
        if (constante.angulo in PESOSUMERGIDO):
            peso_sumergido_x = PESOSUMERGIDO[constante.angulo]
        else:
            peso_sumergido_x = PesoSumergidoX(const, constante.taus, constante.angulo)
            PESOSUMERGIDO[constante.angulo] = peso_sumergido_x
        masa_virtual = MasaVirtual(const, pz_antiguo, vrZ)

        Fx=[fuerza_arrastre_x, peso_sumergido_x, masa_virtual]

        # Fuerzas Y
        fuerza_arrastre_y = Drag(vrY, mgR, coeficiente, const)

        Fy = [fuerza_arrastre_y]

        # Fuerzas Z
        fuerza_arrastre_z = Drag(vrZ, mgR, coeficiente, const)
        peso_sumergido_z = PesoSumergidoZ(const, constante.taus, constante.angulo)
        v_top = VSuperior(vx_antiguo, constante.taus, pz_antiguo, vrY, vrZ)
        v_bot = VBotton(vx_antiguo, constante.taus, pz_antiguo, vrY, vrZ)
        Lift = FuerzaElevacion(const, constante.CL, v_top, v_bot)

        Fz = [fuerza_arrastre_z, peso_sumergido_z, Lift]


        #Nuevas
        vx_actual = Velocidad(vx_antiguo, delta_t, Fx)
        vy_actual = Velocidad(vy_antiguo, delta_t, Fy)
        vz_actual = Velocidad(vz_antiguo, delta_t, Fz)

        px_actual = Posicion(px_antiguo, vx_actual, delta_t)
        py_actual = Posicion(py_antiguo, vy_actual, delta_t)
        pz_actual = Posicion(pz_antiguo, vz_actual, delta_t)

        if(pz_actual < 0.5):
            pz_actual = 0.5

        #nuevas velocidades

        if (pz_actual < pz_antiguo) and (vz_actual*vz_antiguo < 0):
            lista_alturas_alcanzadas.append(pz_antiguo)

        if pz_actual < 0.501:
            vz_actual = vz_actual * -1
            angulo_rebote = AnguloRebote(vz_actual, vx_actual)
            vx_actual = UPrima(angulo_rebote, vz_actual)
            vy_actual = VPrima(vx_actual)
            pz_actual = 0.501
            total_saltos += 1


        nueva_particula.px = px_actual
        nueva_particula.py = py_actual
        nueva_particula.pz = pz_actual
        nueva_particula.vx = vx_actual
        nueva_particula.vy = vy_actual
        nueva_particula.vz = vz_actual

        tiempo_transcurrido += delta_t

        # Fin while

    #print(lista_alturas_alcanzadas)
    h_max = max(lista_alturas_alcanzadas)
    #h_max=1
    h_promedio = numpy.mean(lista_alturas_alcanzadas)
    # Resultado resumido
    resultado = {
        "id":p_id,
        "posiciones_finales":[round(nueva_particula.px,2), round(nueva_particula.py,2), round(nueva_particula.pz,2)],
        "cantidad_de_saltos":total_saltos,
        "altura_maxima":round(h_max,2),
        "promedio_altura_saltos":round(h_promedio,2)
    }
    print(f"Particula #{p_id} simulada en {time.time() - t_lap:.2f} segundos")
    resultados_reales[p_id] = resultado
    return resultado

def GuardarResultadosEnArchivo(nombre_archivo, lista_resultados):
    file = open(f'{nombre_archivo}.out', 'w')
    for r in range(len(lista_resultados)):
        resultado = lista_resultados[r]
        linea = f"{resultado['posiciones_finales'][0]} "
        linea += f"{resultado['posiciones_finales'][1]} "
        linea += f"{resultado['posiciones_finales'][2]} "
        linea += f"{resultado['cantidad_de_saltos']} "
        linea += f"{resultado['altura_maxima']} "
        linea += f"{resultado['promedio_altura_saltos']}"
        linea += "\n"
        file.write(linea)
    file.close()

filename = "input01"
lista_de_particulas = []
datos = ObtenerDatos(f"{filename}.in")
lista_de_particulas = datos[1]

tiempo_total_simulacion = datos[0].tiempo_simulacion
delta_t = datos[0].delta_t
dato = datos[0]


# Valores repetidos
TAUS_TIME_73 = numpy.sqrt(dato.taus) * 73

#Simulacion
#CargarDatos()
const = 1 + dato.relacion_densidad_agua + 0.5
resultados_reales = {}
i = 0
t_cero = time.time()

ejecutar_en_threads = False

if not ejecutar_en_threads:
    # Inicio simulacion en serie
    for p in range(len(lista_de_particulas)): #Simular solo una particula
        particula = lista_de_particulas[p]
        t_lap = time.time()
        SimularParticula(p, particula, dato, const, t_cero)
else:

    # Fin simulacion en serie
    # Inicio simulacion en threads

    threads = list()
    for p in range(len(lista_de_particulas)): #Simular solo una particula
        i += 1
        particula = lista_de_particulas[p]
        t_lap = time.time()
        x = threading.Thread(target=(SimularParticula), args=(p, particula, dato, const, t_cero))
        threads.append(x)
        x.start()

    while(threading.activeCount() > 1): # Esperar que terminen los threads
        pass

    # Fin de simulacion en threads
print(f"{len(lista_de_particulas)} particulas simuladas en {time.time() - t_cero:.2f} segundos")

GuardarResultadosEnArchivo(filename, resultados_reales)
