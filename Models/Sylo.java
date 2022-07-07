package Models;

import java.util.ArrayList;
import java.util.List;

public class Sylo {
    
    // Variables
    double len, wid, dep, d, dt, kn, kt, ktm, floorX, floorZ;
    List<Particle> particles, prevParticles, borders;

    // Constructor
    public Sylo(double len, double wid, double dep, double d, double dt, double ktm) {
        this.len = len;
        this.wid = wid;
        this.dep = dep;
        this.d = d;
        this.dt = dt;
        this.kn = Math.pow(10, 5);
        this.kt = ktm * kn;
        this.floorX = (wid-d)/2;
        this.floorZ = (dep-d)/2;
        this.particles = new ArrayList<>();
        this.prevParticles = new ArrayList<>();
        this.borders = generateBorderParticles();
    }

    public List<Particle> generateBorderParticles() {
        List<Particle> list = new ArrayList<>();
        list.add(new Particle(0, 0, 0, 0, 0, 0, 0.0001, 0));
        list.add(new Particle(wid, 0, 0, 0, 0, 0, 0.0001, 0));
        list.add(new Particle(0, 0, dep, 0, 0, 0, 0.0001, 0));
        list.add(new Particle(wid, 0, dep, 0, 0, 0, 0.0001, 0));
        list.add(new Particle(0, len, 0, 0, 0, 0, 0.0001, 0));
        list.add(new Particle(wid, len, 0, 0, 0, 0, 0.0001, 0));
        list.add(new Particle(0, len, dep, 0, 0, 0, 0.0001, 0));
        list.add(new Particle(wid, len, dep, 0, 0, 0, 0.0001, 0));
        return list;
    }

    public Particle nextEuler(int index) {

        Particle p = particles.get(index);

        Triple<Double, Double, Double> force = Forces.calculateForce(p, this, true, index);

        // Calculamos la posicion en (t - Dt)

        double newVx = p.vx - (dt / p.m) * force.first;
        double newVy = p.vy - (dt / p.m) * force.second;
        double newVz = p.vz - (dt / p.m) * force.third;

        // Calculamos la velocidad en (t - dt)

        double newx = p.x - dt * newVx + (dt*dt) * force.first / (2*p.m);
        double newy = p.y - dt * newVy + (dt*dt) * force.second / (2*p.m);
        double newz = p.z - dt * newVz + (dt*dt) * force.third / (2*p.m);

        return new Particle(newx, newy, newz, newVx, newVy, newVz, p.r, p.m);

    }

    public Particle nextBeeman(int index) {

        // Nos traemos las particulas ...

        Particle current = particles.get(index);
        Particle prev = prevParticles.get(index);

        // Calculamos las fuerzas ...

        Triple<Double, Double, Double> currentForce = Forces.calculateForce(current, this, true, index);
        Triple<Double, Double, Double> prevForce = Forces.calculateForce(prev, this, false, index);

        // Datos a calcular ...

        double rx = 0, ry = 0, rz = 0, newVx = 0, newVy = 0, newVz;

        // Condiciones iniciales

        double r0x = current.x, r0y = current.y, r0z = current.z, v0x = current.vx, v0y = current.vy, v0z = current.vz;

        // Algoritmo de beeman para un paso ...

        // Posiciones
        rx = r0x + (v0x * dt) + (2.0/3.0) * (currentForce.first/current.m) * (dt * dt) - (1.0/6.0) * (prevForce.first/current.m) * (dt * dt);
        ry = r0y + (v0y * dt) + (2.0/3.0) * (currentForce.second/current.m) * (dt * dt) - (1.0/6.0) * (prevForce.second/current.m) * (dt * dt);
        rz = r0z + (v0z * dt) + (2.0/3.0) * (currentForce.third/current.m) * (dt * dt) - (1.0/6.0) * (prevForce.third/current.m) * (dt * dt);


        // Predicciones
        double predVx = v0x + (3.0/2.0) * (currentForce.first/current.m) * dt - 0.5 * (prevForce.first/current.m) * dt;
        double predVy = v0y + (3.0/2.0) * (currentForce.second/current.m) * dt - 0.5 * (prevForce.second/current.m) * dt;
        double predVz = v0z + (3.0/2.0) * (currentForce.third/current.m) * dt - 0.5 * (prevForce.third/current.m) * dt;


        // Correcciones
        Particle newp = new Particle(rx, ry, rz, predVx, predVy, predVz, current.r, current.m);
        Triple<Double, Double, Double> newForce = Forces.calculateForce(newp, this,true, index);

        newVx= v0x + (1.0/3.0) * (newForce.first/current.m) * dt + (5.0/6.0) * (currentForce.first/current.m) * dt - (1.0/6.0) * (prevForce.first/current.m) * dt;
        newVy= v0y + (1.0/3.0) * (newForce.second/current.m) * dt + (5.0/6.0) * (currentForce.second/current.m) * dt - (1.0/6.0) * (prevForce.second/current.m) * dt; 
        newVz= v0z + (1.0/3.0) * (newForce.third/current.m) * dt + (5.0/6.0) * (currentForce.third/current.m) * dt - (1.0/6.0) * (prevForce.third/current.m) * dt; 

        return new Particle(rx, ry, rz, newVx, newVy, newVz, current.r, current.m);
    }


}
