package Models;

import javax.swing.text.html.HTMLDocument.RunElement;

public class Particle {
    
    // Variables
    double x, y, z, vx, vy, vz, r, m;
    boolean ignore;

    // Constructor
    public Particle (double x, double y, double z, double vx, double vy, double vz, double r, double m) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.vx = vx;
        this.vy = vy;
        this.vz = vz;
        this.r = r;
        this.m = m;
        this.ignore = false;
    }

    public double getVelocity() { return Math.sqrt(vx*vx + vy*vy + vz*vz); }
    public double getDistance(Particle o) { return Math.sqrt( (this.x - o.x)*(this.x - o.x) + (this.y - o.y)*(this.y - o.y) + (this.z - o.z)*(this.z - o.z) ); }
    public double getOverlap(Particle o) { return this.r + o.r - getDistance(o); }
    public Triple<Double, Double, Double> getRV(Particle o) {
        double rvx = this.vx - o.vx;
        double rvy = this.vy - o.vy;
        double rvz = this.vz - o.vz;
        Triple<Double, Double, Double> xyz = new Triple<Double, Double, Double>(rvx, rvy, rvz);
        return xyz;
    }

}
