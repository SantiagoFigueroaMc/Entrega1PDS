import math
import random
from random import seed
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import time

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
    resultado = -0.75 * (1/var) * coeficiente_de_arrastre * v_relativa * magnitud_v_relativa
    return resultado


def VRelativaX(v_antigua, taus,pz):
    usando = VFluido(taus,pz)
    return v_antigua - usando

def VFluido(taus,pz):
    if 73 * math.sqrt(taus) < 5:
        usando = 2.5 * math.log(73 * math.sqrt(taus) * abs(pz)) + 5.5
    elif 5 <= 73 * math.sqrt(taus) < 70:
        usando = 2.5 * math.log(73 * math.sqrt(taus) * abs(pz)) + 5.5 - 2.5 * math.log( 1 + 0.3 * 73 * math.sqrt(taus))
    else:
        usando = 2.5 * math.log(30 * abs(pz))

    return usando

def MagnitudVRelativa(u_relativa, v_relativa, w_relativa):
    resultado = math.sqrt( (u_relativa**2) + (v_relativa**2) + (w_relativa**2) )
    #print(f"MagnitudVRelativa: {resultado:.2f} (u:{u_relativa:.2f},v:{v_relativa:.2f},w:{w_relativa:.2f})")
    # TODO: Sacar condicion de minimo
    #return min(10, resultado)
    return resultado

def CoeficienteDeArrastre(rep):
    componente1 = 1 + (10**4) * (rep**-0.5)
    division1 = 0.208 / componente1
    componente2 = rep * ( 1 + 0.15 * math.sqrt(rep) + 0.017 * rep)
    resultado = 24 / (componente2 - division1)
    return resultado

def Rep(magnitud_v_relativa, taus):
    resultado = magnitud_v_relativa * math.sqrt(taus) * 73
    return resultado

def PesoSumergidoX(constante, taus, teta):
    resultado = (1/constante) * math.sin(teta) * (1/taus)
    return resultado

def PesoSumergidoZ(constante, taus, teta):
    resultado = - (1/constante) * math.cos(teta) * (1/taus)
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
    resultado = math.atan(w_prima/u_actual)
    while resultado > math.radians(75):
        resultado = math.atan(w_prima/u_actual)
    return resultado

def UPrima(alfa,w_prima):
    error = random.uniform(math.radians(0), math.radians(10))
    resultado = w_prima / math.tan(alfa + error)
    return resultado

def VPrima(u_prima):
    alfa = random.uniform(math.radians(-10), math.radians(10))
    resultado = u_prima * math.tan(alfa)
    return resultado

def Posicion(p_vieja, v_actual, delta_t):
    resultado = p_vieja + v_actual * delta_t
    return resultado

def Velocidad(v_vieja, delta_t, fuerzas):
    sumatoria = 0
    for fuerza in fuerzas:
        sumatoria += fuerza
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

def SimularParticula(particula_simulada,constante,const):

    # Para plot
    global t
    t = []

    global lista_posiciones_x
    global lista_posiciones_y
    global lista_posiciones_z
    lista_posiciones_x = []
    lista_posiciones_y = []
    lista_posiciones_z = []

    global lista_velocidades_x
    global lista_velocidades_y
    global lista_velocidades_z
    lista_velocidades_x = []
    lista_velocidades_y = []
    lista_velocidades_z = []

    global lista_suma_fuerzas_x
    global lista_suma_fuerzas_y
    global lista_suma_fuerzas_z
    lista_suma_fuerzas_x = []
    lista_suma_fuerzas_y = []
    lista_suma_fuerzas_z = []

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


        # TODO: Simular
        px_antiguo = nueva_particula.px
        py_antiguo = nueva_particula.py
        pz_antiguo = nueva_particula.pz
        vx_antiguo = nueva_particula.vx
        vy_antiguo = nueva_particula.vy
        vz_antiguo = nueva_particula.vz

        # General
        vrX = VRelativaX(vx_antiguo, constante.taus,pz_antiguo)
        vrY = vy_antiguo
        vrZ = vz_antiguo

        mgR = MagnitudVRelativa(vrX, vrY, vrZ)
        rep = Rep(mgR, constante.taus)
        coeficiente =  CoeficienteDeArrastre(rep)

        # Fuerzas X
        fuerza_arrastre_x = Drag(vrX,mgR,coeficiente,const)


        peso_sumergido_x = PesoSumergidoX(const, constante.taus, constante.angulo)

        masa_virtual = MasaVirtual(const, pz_antiguo, vrZ)

        Fx=[fuerza_arrastre_x, peso_sumergido_x, masa_virtual]

        # Fuerzas Y
        fuerza_arrastre_y = Drag(vrY,mgR,coeficiente,const)
        Fy = [fuerza_arrastre_y]



        # Fuerzas Z
        fuerza_arrastre_z = Drag(vrZ,mgR,coeficiente,const)

        peso_sumergido_z = PesoSumergidoZ(const, constante.taus, constante.angulo)

        v_top = VSuperior(vx_antiguo,constante.taus,pz_antiguo, vrY, vrZ)

        v_bot = VBotton(vx_antiguo,constante.taus,pz_antiguo, vrY, vrZ)


        Lift = FuerzaElevacion(const,constante.CL, v_top, v_bot)

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

        #direccion_z_antiguo = direccion_z
        #direccion_z = pz_actual > pz_antiguo

        #print(direccion_z, direccion_z_antiguo)
        #print(pz_actual, pz_antiguo)
        datos = f"Fuerzas:\n"
        datos+= f"    Fx: (Arrastre:{Fx[0]:.2f}, Sumergido:{Fx[1]:.2f}, MasaVirtual:{Fx[2]:.2f})\n"
        datos+= f"    Fy: (Arrastre:{Fy[0]:.2f})\n"
        datos+= f"    Fz: (Arrastre:{Fz[0]:.2f}, Sumergido:{Fz[1]:.2f}, Lift:{Fz[2]:.2f})\n"
        datos+= f"Velocidades:\n"
        datos+= f"    Vx: ({vx_actual:.2f})\n"
        datos+= f"    Vy: ({vy_actual:.2f})\n"
        datos+= f"    Vz: ({vz_actual:.2f})\n"
        datos+= f"Pocisiones:\n"
        datos+= f"    Px: ({px_actual:.2f})\n"
        datos+= f"    Py: ({py_actual:.2f})\n"
        datos+= f"    Pz: ({pz_actual:.2f})\n"
        #if (pz_actual < 0.501):
        #print(datos)

        if (pz_actual < pz_antiguo) and (vz_actual*vz_antiguo < 0):
            #print(pz_antiguo)
            lista_alturas_alcanzadas.append(pz_antiguo)

        if pz_actual < 0.501:
            #print("hola")
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


        #Sacar para optimizar
        #para pyplot
        lista_posiciones_x.append(nueva_particula.px)
        lista_posiciones_y.append(nueva_particula.py)
        lista_posiciones_z.append(nueva_particula.pz)

        lista_velocidades_x.append(nueva_particula.vx)
        lista_velocidades_y.append(nueva_particula.vy)
        lista_velocidades_z.append(nueva_particula.vz)

        lista_suma_fuerzas_x.append(sum(Fx))
        lista_suma_fuerzas_y.append(sum(Fy))
        lista_suma_fuerzas_z.append(sum(Fz))
        t.append(tiempo_transcurrido)

        tiempo_transcurrido += delta_t

    #print(lista_alturas_alcanzadas)
    h_max = max(lista_alturas_alcanzadas)
    #h_max=1
    h_promedio = numpy.mean(lista_alturas_alcanzadas)
    # Resultado resumido
    resultado = {
        "posiciones_finales":[round(nueva_particula.px,2), round(nueva_particula.py,2), round(nueva_particula.pz,2)],
        "cantidad_de_saltos":total_saltos,
        "altura_maxima":round(h_max,2),
        "promedio_altura_saltos":round(h_promedio,2)
    }
    return resultado

def GuardarResultadosEnArchivo(nombre_archivo, lista_resultados):
    file = open(f'{nombre_archivo}.out', 'w')
    for resultado in lista_resultados:
        linea = f"{resultado['posiciones_finales'][0]} "
        linea += f"{resultado['posiciones_finales'][1]} "
        linea += f"{resultado['posiciones_finales'][2]} "
        linea += f"{resultado['cantidad_de_saltos']} "
        linea += f"{resultado['altura_maxima']} "
        linea += f"{resultado['promedio_altura_saltos']}"
        linea += "\n"
        file.write(linea)
    file.close()

filename = "input02"
lista_de_particulas = []
datos = ObtenerDatos(f"{filename}.in")
lista_de_particulas = datos[1]

tiempo_total_simulacion = datos[0].tiempo_simulacion
delta_t = datos[0].delta_t
dato = datos[0]

#Simulacion
#CargarDatos()
const = 1 + dato.relacion_densidad_agua + 0.5
resultados_reales =[]
i = 0
t_cero = time.time()
for particula in lista_de_particulas: #Simular solo una particula
      i += 1
      t_lap = time.time()
      resultados_reales.append(SimularParticula(particula,dato,const))
      print(f"particula #{i} in {time.time() - t_lap:.2f} seconds ({time.time() - t_cero:.2f} s total)")
      #print(resultados_reales)

t = numpy.arange(0.0, len(lista_posiciones_x)*delta_t, delta_t)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)

ax1.plot(t, lista_posiciones_x)
ax1.set(xlabel='t', ylabel='X', title='pX in time')
ax1.grid()

ax2.plot(t, lista_posiciones_y)
ax2.set(xlabel='t', ylabel='Y', title='pY in time')
ax2.grid()

promedio_altura_saltos_lista = numpy.empty(len(t))
promedio_altura_saltos_lista.fill(resultados_reales[0]['promedio_altura_saltos'])
ax3.plot(t, lista_posiciones_z, t, promedio_altura_saltos_lista)
ax3.set(xlabel='t', ylabel='Z', title='pZ in time')
ax3.grid()

ax4.plot(t, lista_velocidades_x)
ax4.set(xlabel='t', ylabel='X', title='vX in time')
ax4.grid()

ax5.plot(t, lista_velocidades_y)
ax5.set(xlabel='t', ylabel='Y', title='vY in time')
ax5.grid()

ax6.plot(t, lista_velocidades_z)
ax6.set(xlabel='t', ylabel='Z', title='vZ in time')
ax6.grid()

ax7.plot(t, lista_suma_fuerzas_x)
ax7.set(xlabel='t', ylabel='X', title='FX in time')
ax7.grid()

ax8.plot(t, lista_suma_fuerzas_y)
ax8.set(xlabel='t', ylabel='Y', title='FY in time')
ax8.grid()

ax9.plot(t, lista_suma_fuerzas_z)
ax9.set(xlabel='t', ylabel='Z', title='FZ in time')
ax9.grid()

#fig.tight_layout()
plt.show()

ax10 = plt.axes(projection = '3d')
z = lista_posiciones_z
y = lista_posiciones_y
x = lista_posiciones_x

ax10.plot3D(x[:5000], y[:5000], z[:5000],color='green')

ax10.set_title('xyz position in time')
plt.show()


GuardarResultadosEnArchivo(filename, resultados_reales)
