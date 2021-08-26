import numpy as np
from Ecuaciones import *
from Archivo import *
from Modelos import *
import time

archivo="c6-p2"

def ObtenerDatos(nombre_archivo):
    contador_linea=0
    pos_x = []
    pos_y = []
    pos_z = []
    vel_x = []
    vel_y = []
    vel_z = []
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
                pos_x.append(float(linea[0]))
                pos_y.append(float(linea[1]))
                pos_z.append(float(linea[2]))
                vel_x.append(float(linea[3]))
                vel_y.append(float(linea[4]))
                vel_z.append(float(linea[5]))

            contador_linea+=1
        else:
            break

    archivo.close()
    return simulacion,np.array(pos_x),np.array(pos_y),np.array(pos_z),np.array(vel_x),np.array(vel_y),np.array(vel_z)

# posiciones y velocidades para tres particulas
entorno,pos_x, pos_y, pos_z, vel_x, vel_y, vel_z = ObtenerDatos(f"{archivo}.in")
t_simulacion = entorno.tiempo_simulacion
delta_t = entorno.delta_t
const = 1 + entorno.relacion_densidad_agua + 0.5
lista_alturas_alcanzadas = []

for l in range(len(pos_x)):
    lista_alturas_alcanzadas.append([])

total_saltos = [0] * len(pos_x)
t_cero = time.time()



# simular 100 iteraciones
while(t_simulacion > 0): 
    #print("Tiempo restante:", t_simulacion)
    vel_antiguo = vel_z
    pos_antiguo = pos_z
    v_rel_x = VRelativaX(vel_x, entorno.taus, pos_z)
    
    mgR =  MagnitudVRelativa(v_rel_x, vel_y, vel_z)
    rep = Rep(mgR, entorno.taus)
    coeficiente = CoeficienteDeArrastre(rep)

    # Fuerzas en X
    arrastre_x = Drag(v_rel_x, mgR, coeficiente, const)
    sumergido_x = PesoSumergidoX(const, entorno.taus, entorno.angulo)
    masa_virtual = MasaVirtual(const, pos_z, vel_z)
    Fx = [arrastre_x, sumergido_x, masa_virtual]

    # Fuerzas en Y
    arrastre_y = Drag(vel_y, mgR, coeficiente, const)
    Fy = [arrastre_y]

    # Fuerzas en Z
    arrastre_z = Drag(vel_z, mgR, coeficiente, const)
    sumergido_z = PesoSumergidoZ(const, entorno.taus, entorno.angulo)
    v_top = VSuperior(vel_x, entorno.taus, pos_z, vel_y, vel_z)
    v_bot = VBotton(vel_x, entorno.taus, pos_z, vel_y, vel_z)
    lift = FuerzaElevacion(const, entorno.CL, v_top, v_bot)
    Fz = [arrastre_z, sumergido_z, lift]

    # Nuevas velocidades
    vel_x = Velocidad(vel_x, delta_t, Fx)
    vel_y = Velocidad(vel_y, delta_t, Fy)
    vel_z = Velocidad(vel_z, delta_t, Fz)
    # Nuevas posiciones
    pos_x = Posicion(pos_x, vel_x, delta_t)
    pos_y = Posicion(pos_y, vel_y, delta_t)
    pos_z = Posicion(pos_z, vel_z, delta_t)
    

    for p in np.where(pos_z < pos_antiguo)[0]:
        lista_alturas_alcanzadas[p].append(pos_antiguo[p])

    for p in np.where(pos_z < 0.5)[0]:
        vel_z[p] = vel_z[p] * -1
        angulo_rebote = AnguloRebote(vel_z[p], vel_x[p])
        vel_x[p] = UPrima(angulo_rebote, vel_z[p])
        vel_y[p] = VPrima(vel_x[p])
        pos_z[p] = 0.501
        total_saltos[p] += 1

    t_simulacion -= delta_t

resultados=[]


for p in range(len(pos_x)):
    h_max = max(lista_alturas_alcanzadas[p])
    #h_max=1
    h_promedio = np.mean(lista_alturas_alcanzadas[p])
    resultado = {
        "posiciones_finales":[round(pos_x[p],2), round(pos_y[p],2), round(pos_z[p],2)],
        "cantidad_de_saltos":total_saltos[p],
        "altura_maxima":round(h_max,2),
        "promedio_altura_saltos":round(h_promedio,2)
    }
    resultados.append(resultado)

#print(resultados)
GuardarResultadosEnArchivo(archivo, resultados)


print(f"Done in {time.time() - t_cero:.2f} seconds")
#print(pos_z)
