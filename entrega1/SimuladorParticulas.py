from Ecuaciones import *
from random import seed
import time
import multiprocessing as mp
import numpy as np
from Archivo import ObtenerDatos, GuardarResultadosEnArchivo

def SimularParticula(p_id, particula_simulada, entorno, resultados_reales):
    print(f"Simulando particula #{p_id}")
    t_lap = time.time()

    # Variables necesarias
    tiempo_transcurrido = 0
    tiempo_total_simulacion = entorno.tiempo_simulacion
    delta_t = entorno.delta_t

    total_saltos = 0
    lista_alturas_alcanzadas = []
    h_max = 0.0
    h_promedio = 0.0

    direccion_z = False #True sube, False baja
    direccion_z_antiguo = False

    nueva_particula = particula_simulada

    # Variables que se calculan repetidas veces
    const = 1 + entorno.relacion_densidad_agua + 0.5

    # Diccionarios para memorizacion (individial)
    REP = {}
    PESOSUMERGIDO = {}


    # Simulacion
    while (tiempo_transcurrido < tiempo_total_simulacion):
        px_antiguo = nueva_particula.px
        py_antiguo = nueva_particula.py
        pz_antiguo = nueva_particula.pz
        vx_antiguo = nueva_particula.vx
        vy_antiguo = nueva_particula.vy
        vz_antiguo = nueva_particula.vz

        # General
        vrX = VRelativaX(vx_antiguo, entorno.taus, pz_antiguo)
        vrY = vy_antiguo
        vrZ = vz_antiguo

        mgR = MagnitudVRelativa(vrX, vrY, vrZ)
        rep = Rep(mgR, entorno.taus)
        if (rep in REP):
            coeficiente = REP[rep]
        else:
            coeficiente =  CoeficienteDeArrastre(rep)
            REP[rep] = coeficiente

        # Fuerzas X
        fuerza_arrastre_x = Drag(vrX, mgR, coeficiente, const)
        if (entorno.angulo in PESOSUMERGIDO):
            peso_sumergido_x = PESOSUMERGIDO[entorno.angulo]
        else:
            peso_sumergido_x = PesoSumergidoX(const, entorno.taus, entorno.angulo)
            PESOSUMERGIDO[entorno.angulo] = peso_sumergido_x
        masa_virtual = MasaVirtual(const, pz_antiguo, vrZ)

        Fx=[fuerza_arrastre_x, peso_sumergido_x, masa_virtual]

        # Fuerzas Y
        fuerza_arrastre_y = Drag(vrY, mgR, coeficiente, const)

        Fy = [fuerza_arrastre_y]

        # Fuerzas Z
        fuerza_arrastre_z = Drag(vrZ, mgR, coeficiente, const)
        peso_sumergido_z = PesoSumergidoZ(const, entorno.taus, entorno.angulo)
        v_top = VSuperior(vx_antiguo, entorno.taus, pz_antiguo, vrY, vrZ)
        v_bot = VBotton(vx_antiguo, entorno.taus, pz_antiguo, vrY, vrZ)
        Lift = FuerzaElevacion(const, entorno.CL, v_top, v_bot)

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
    h_promedio = np.mean(lista_alturas_alcanzadas)
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

def main():
    manager = mp.Manager()
    results = manager.dict()

    filename = "input03"

    total_time = time.time()

    entorno, particulas = ObtenerDatos(f"{filename}.in")
    
    particulas = particulas[:60] # Pool no aguanta más de 62 trabajadores

    pool = mp.Pool(len(particulas))

    jobs = []

    for p in range(len(particulas)):
        particula = particulas[p]
        job = pool.apply_async(SimularParticula, (p, particula, entorno, results))
        jobs.append(job)

    for job in jobs:
        job.get()

    pool.close()
    pool.join()

    lista_resultados = []
    for i in range(len(particulas)):
        lista_resultados.append(results[i])

    print(f"{len(particulas)} particula(s) simladas en {time.time() - total_time:.2f} segundos")
    GuardarResultadosEnArchivo(filename, lista_resultados)

if __name__ == "__main__":
    main()
