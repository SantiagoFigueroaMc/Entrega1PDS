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
