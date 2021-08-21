import numpy as np
from Ecuaciones import *
from Archivo import *
from Modelos import *
import time

entorno = Variables
entorno.taus = 0.067
entorno.angulo = 0.001
entorno.relacion_densidad_agua = 1.65
entorno.CL = 0.2

t_simulacion = 100
delta_t = 0.001

const = 1 + entorno.relacion_densidad_agua + 0.5
    
# posiciones para tres particulas
pos_x = np.array([0, 0, 0, 0, 0])
pos_y = np.array([0, 1, 2, 3, 4])
pos_z = np.array([.6, .6, .6, .6, .6])
# velocidades para tres particulas
vel_x = np.array([4.15 , 3.15, 2.15, 3.15, 4.15])
vel_y = np.array([0.1, -0.1, 0.1, -0.1, 0.1])
vel_z = np.array([1.85, 2, 2, 1.5, 1.2])

t_cero = time.time()
# simular 100 iteraciones
while(t_simulacion > 0): 
    #print("Tiempo restante:", t_simulacion)
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

    for p in range(len(pos_z)):
        if pos_z[p] < 0.5:
            #print("Rebote", pos_z[p])
            pos_z[p] = 0.501
            vel_z[p] = vel_z[p] * -1
            angulo_rebote = AnguloRebote(vel_z[p], vel_x[p])
            vel_x[p] = UPrima(angulo_rebote, vel_z[p])
            vel_y[p] = VPrima(vel_x[p])
            #total_saltos += 1

    t_simulacion -= delta_t

print(f"Done in {time.time() - t_cero:.2f} sconds")
print(pos_z)
