from Modelos import Particula, Variables

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
