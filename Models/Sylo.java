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
        list.add(new Particle(0, 0, 0, 0, 0, 0, 0.0001, 0, false));
        list.add(new Particle(wid, 0, 0, 0, 0, 0, 0.0001, 0, false));
        list.add(new Particle(0, 0, dep, 0, 0, 0, 0.0001, 0, false));
        list.add(new Particle(wid, 0, dep, 0, 0, 0, 0.0001, 0, false));
        list.add(new Particle(0, len, 0, 0, 0, 0, 0.0001, 0, false));
        list.add(new Particle(wid, len, 0, 0, 0, 0, 0.0001, 0, false));
        list.add(new Particle(0, len, dep, 0, 0, 0, 0.0001, 0, false));
        list.add(new Particle(wid, len, dep, 0, 0, 0, 0.0001, 0, false));
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

        return new Particle(newx, newy, newz, newVx, newVy, newVz, p.r, p.m, false);

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
        Particle newp = new Particle(rx, ry, rz, predVx, predVy, predVz, current.r, current.m, current.ignore);
        Triple<Double, Double, Double> newForce = Forces.calculateForce(newp, this,true, index);

        newVx= v0x + (1.0/3.0) * (newForce.first/current.m) * dt + (5.0/6.0) * (currentForce.first/current.m) * dt - (1.0/6.0) * (prevForce.first/current.m) * dt;
        newVy= v0y + (1.0/3.0) * (newForce.second/current.m) * dt + (5.0/6.0) * (currentForce.second/current.m) * dt - (1.0/6.0) * (prevForce.second/current.m) * dt; 
        newVz= v0z + (1.0/3.0) * (newForce.third/current.m) * dt + (5.0/6.0) * (currentForce.third/current.m) * dt - (1.0/6.0) * (prevForce.third/current.m) * dt; 

        return new Particle(rx, ry, rz, newVx, newVy, newVz, current.r, current.m, current.ignore);
    }

    public Particle placeNewParticle(int seconds) {
        long currentTime= System.currentTimeMillis();
        long end = currentTime + (seconds * 1000);

        double radiusLow = 0.04/2;
        double radiusHigh = 0.05/2;

        double y_high = len;
        double y_low = (4 * len)/5;
        double x_high = wid;
        double z_high = dep;
        
        while((currentTime = System.currentTimeMillis()) < end) {

            double rand_r = (Math.random() * (radiusHigh-radiusLow)) + radiusLow;
            double rand_x = Math.random() * x_high;
            double rand_z = Math.random() * z_high;
            double rand_y = Math.random() * (y_high - y_low) + y_low;

            Particle p = new Particle(rand_x, rand_y, rand_z, 0, 0, 0, rand_r, 0.01, false);

            while(true) {
                Particle sup = new Particle(p.x, len, p.z, 0, 0, 0, 0, 0, false);
                Particle der = new Particle(wid, p.y, p.z, 0, 0, 0, 0, 0, false);
                Particle izq = new Particle(0, p.y, p.z, 0, 0, 0, 0, 0, false);
                Particle front = new Particle(p.x, p.y, dep, 0, 0, 0, 0, 0, false);
                Particle back = new Particle(p.x, p.y, 0, 0, 0, 0, 0, 0, false);
                Particle inf = new Particle(p.x, 0, p.z, 0, 0, 0, 0, 0, false);

                if(p.getOverlap(sup) >= 0 || p.getOverlap(inf) >= 0 || p.getOverlap(izq) >= 0 || p.getOverlap(der) >= 0 || p.getOverlap(front) >= 0 || p.getOverlap(back) >= 0){
                    rand_x = Math.random() * x_high;
                    rand_z = Math.random() * z_high;
                    rand_y = Math.random() * (y_high - y_low) + y_low;
                    p = new Particle(rand_x, rand_y, rand_z, 0, 0, 0, rand_r, p.m, false);
                } else {
                    break;
                }
            }

            boolean flag = true;
            for(Particle other : particles){
                if(p.getOverlap(other) >= 0){
                    flag = false;
                }
            }
            if(flag) {
                return p;
            }
        }
        return null;
    }

    public void simulateUniverse(int seconds) {

        List<Boolean> first = new ArrayList<>();

        for(Particle p : particles)
            first.add(true);

        int step = 1;
        int n = 50000;
        int tOutput = 1000;

        while((dt*step) < 5) {
            int index = 0;
            for(Particle p : particles) {
                if(first.get(index)) {
                    prevParticles.set(index, nextEuler(index));
                    first.set(index, false);
                }
                Particle aux = particles.get(index);
                Particle newParticle = nextBeeman(index);
                //chequear si la particula se fue del sylo la vuelvo a poner
                if(newParticle.x > wid || newParticle.x < 0 || newParticle.y > len || newParticle.z < 0 || newParticle.z > dep){
                    first.set(index, true);
                    Particle ret = placeNewParticle(seconds);
                    particles.set(index, ret);
                    // reincerciones no deseadas
                }
                else if(newParticle.y <= -len/10) {
                    first.set(index, true);
                    Particle ret = placeNewParticle(seconds);
                    particles.set(index, ret);
                } else {
                    particles.set(index, newParticle);
                    prevParticles.set(index, aux);
                }
                index ++;
                if(step % tOutput == 0) {
                    OutputParser.writeUniverse(particles, borders, step*dt);
                }
            }
            if((step) % n ==  0)
                System.out.println(step*dt);
            step++;
        }
    }

    public void simulateEJ1_2(int seconds, String fn) {

        List<Boolean> first = new ArrayList<>();

        for(Particle p : particles)
            first.add(true);

        int step = 1;
        int n = 50000;
        int printcsv = 10000;
        int escape_counter = 0;
        while((dt*step) < 5) {
            int index = 0;
            for(Particle p : particles) {
                if(first.get(index)) {
                    prevParticles.set(index, nextEuler(index));
                    first.set(index, false);
                }
                Particle aux = particles.get(index);
                Particle newParticle = nextBeeman(index);
                if(newParticle.x > wid || newParticle.x < 0 || newParticle.y > len || newParticle.z < 0 || newParticle.z > dep){
                    escape_counter++;
                    first.set(index, true);
                    Particle ret = placeNewParticle(seconds);
                    particles.set(index, ret);
                    // reincerciones no deseadas
                }
                else if(newParticle.y <= -len/10) {
                    escape_counter++;
                    first.set(index, true);
                    Particle ret = placeNewParticle(seconds);
                    particles.set(index, ret);
                } else {
                    particles.set(index, newParticle);
                    prevParticles.set(index, aux);
                }
                index ++;
            }
            if(step % printcsv == 0){
                OutputParser.parseEj1_2(step*dt, escape_counter, fn);
                escape_counter = 0;
            }
            if((step) % n ==  0)
                System.out.println(step*dt);
            step++;
        }
    }

    public void simulateEJ3_4(int seconds, String fn) {

        List<Boolean> first = new ArrayList<>();

        for(Particle p : particles)
            first.add(true);

        int step = 1;
        int n = 50000;
        int printcsv = 10000;
        double ke = 0;
        OutputParser.parseEj3_4(0, 0, fn);
        while((dt*step) < 10) {
            int index = 0;
            for(Particle p : particles) {
                if(first.get(index)) {
                    prevParticles.set(index, nextEuler(index));
                    first.set(index, false);
                }
                Particle aux = particles.get(index);
                Particle newParticle = nextBeeman(index);
                if(newParticle.x > wid || newParticle.x < 0 || newParticle.y > len || newParticle.z < 0 || newParticle.z > dep){
                    first.set(index, true);
                    Particle ret = placeNewParticle(seconds);
                    particles.set(index, ret);
                    // reincerciones no deseadas
                }
                else if(newParticle.y <= -len/10) {
                    first.set(index, true);
                    Particle ret = placeNewParticle(seconds);
                    particles.set(index, ret);
                } else {
                    ke += 0.5 * newParticle.m * Math.pow(newParticle.getVelocity(),2);
                    particles.set(index, newParticle);
                    prevParticles.set(index, aux);
                }
                index ++;
            }
            if(step % printcsv == 0){
                OutputParser.parseEj3_4(step*dt, ke, fn);
            }
            ke = 0;
            if((step) % n ==  0)
                System.out.println(step*dt);
            step++;
        }
    }

    public void populate(double seconds) {

        long currentTime= System.currentTimeMillis();
        long end = currentTime + (int)(seconds * 1000);

        double radiusLow = 0.04/2;
        double radiusHigh = 0.05/2;

        double y_high = len;
        double x_high = wid;
        double z_high = dep;
        boolean first = true;
        int i = 0;
        // while((currentTime = System.currentTimeMillis()) < end) {
        while(i < 300) {

            double rand_r = (Math.random() * (radiusHigh-radiusLow)) + radiusLow;
            double rand_x = Math.random() * x_high;
            double rand_y = Math.random() * y_high;
            double rand_z = Math.random() * z_high;

            Particle p = new Particle(rand_x, rand_y, rand_z, 0, 0, 0, rand_r, 0.01, false);

            while(true) {
                Particle sup = new Particle(p.x, len, p.z, 0, 0, 0, 0, 0, false);
                Particle der = new Particle(wid, p.y, p.z, 0, 0, 0, 0, 0, false);
                Particle izq = new Particle(0, p.y, p.z, 0, 0, 0, 0, 0, false);
                Particle front = new Particle(p.x, p.y, dep, 0, 0, 0, 0, 0, false);
                Particle back = new Particle(p.x, p.y, 0, 0, 0, 0, 0, 0, false);
                Particle inf = new Particle(p.x, 0, p.z, 0, 0, 0, 0, 0, false);

                if(p.getOverlap(sup) >= 0 || p.getOverlap(inf) >= 0 || p.getOverlap(izq) >= 0 || p.getOverlap(der) >= 0 || p.getOverlap(front) >= 0 || p.getOverlap(back) >= 0){
                    rand_x = Math.random() * x_high;
                    rand_y = Math.random() * y_high;
                    rand_z = Math.random() * z_high;
                    p = new Particle(rand_x, rand_y, rand_z, 0, 0, 0, rand_r, p.m, false);
                } else {
                    break;
                }
            }
            if(first) {
                particles.add(p);
                i++;
                first = false;
            } else {
                boolean flag = true;
                for(Particle other : particles){
                    if(p.getOverlap(other) >= 0){
                        flag = false;
                    }
                }
                if(flag){
                    particles.add(p);
                    i++;
                }
            }
        }

        int in = 0;
        for(Particle part : particles){
            prevParticles.add(nextEuler(in));
            in++;
        }
        
        OutputParser.writeUniverse(particles, borders, 0);
        System.out.println("Finalizo el populate - " + particles.size() + " Particulas");
    }
}
