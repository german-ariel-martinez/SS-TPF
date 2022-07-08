package Models;

import java.util.List;

public class Forces {

    // Metodos
    public static Triple<Double, Double, Double> getPerpendicularVersor(Triple<Double, Double, Double> v1, Triple<Double, Double, Double> v2){
        double nx = v1.second * v2.third - v1.third * v2.second;
        double ny = v2.first * v1.third - v2.third * v1.first;
        double nz = v1.first * v2.second - v1.second * v2.first;
        return new Triple<Double,Double,Double>(nx, ny, nz);
    }

    public static double getFN(Particle p1, Particle p2, Sylo s) {
        return - s.kn * p1.getOverlap(p2);
    }

    public static double getFT(Particle p1, Particle p2, Sylo s) {

        double xDiff = p2.x - p1.x;
        double yDiff = p2.y - p1.y;
        double zDiff = p2.z - p1.z;
        double distance = p1.getDistance(p2);
        double enx = xDiff/distance;
        double eny = yDiff/distance;
        double enz = zDiff/distance;

        // Calculamos los versores tangenciales y normales ...
        Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
        Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);

        // Nos traemos la velocidad relativa entre las particulas ...
        Triple<Double, Double, Double> v = p1.getRV(p2);

        // Solo nos interesa la componente tangencial ...
        double term = ( v.first * et.first + v.second * et.second + v.third * et.third);

        // Retornamos
        return - s.kt * p1.getOverlap(p2) * term;

    }

    public static double getFT2(Particle p1, Particle p2, Sylo s) {

        double xDiff = p2.x - p1.x;
        double yDiff = p2.y - p1.y;
        double zDiff = p2.z - p1.z;
        double distance = p1.getDistance(p2);
        double enx = xDiff/distance;
        double eny = yDiff/distance;
        double enz = zDiff/distance;

        // Calculamos los versores tangenciales y normales ...
        Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
        Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
        Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);

        // Nos traemos la velocidad relativa entre las particulas ...
        Triple<Double, Double, Double> v = p1.getRV(p2);

        // Solo nos interesa la componente tangencial ...
        double term = ( v.first * et2.first + v.second * et2.second + v.third * et2.third);

        // Retornamos
        return - s.kt * p1.getOverlap(p2) * term;

    }

    public static Triple<Double, Double, Double> calculateForce(Particle p, Sylo s, boolean current, int index) {

        // Las componentes de la fuerza ...
        double fxT = 0, fyT = 0, fzT = 0;
        double fn = 0, ft = 0, ft2 = 0;

        // En que array trabajamos

        List<Particle> array;

        if(current)
            array = s.particles;
        else
            array = s.prevParticles;

        // Comparo contra las particulas ...
        for (int i = 0; i < array.size(); i++) {
            if(i != index){
                if(p.getOverlap(array.get(i)) > 0) {
                    // Nos traemos FN y FT
                    fn = getFN(p, array.get(i), s);
                    ft = getFT(p, array.get(i), s);
                    ft2 = getFT2(p, array.get(i), s);
                    // Calculamos las fuerzas totales en X e Y ...
                    double xDiff = array.get(i).x - p.x;
                    double yDiff = array.get(i).y - p.y;
                    double zDiff = array.get(i).z - p.z;
                    double distance = p.getDistance(array.get(i));
                    double enx = xDiff/distance;
                    double eny = yDiff/distance;
                    double enz = zDiff/distance;
                    Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
                    Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
                    Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);
                    fxT += fn * enx + ft * (-eny) + ft2 * et2.first;
                    fyT += fn * eny + ft * enx + ft2 * et2.second;
                    fzT += fn * enz + ft * enz + ft2 * et2.third;
                    if(fxT > 1000 || fyT > 1000 || fzT > 1000){
                        System.out.println("ME fui a la mierda entre parts");
                        System.exit(1);
                    }
                }
            }
        }

        // Creamos las paredes que en este caso son particulas ...
        Particle sup = new Particle(p.x, s.len, p.z, 0, 0, 0, 0, 0);
        Particle der = new Particle(s.wid, p.y, p.z, 0, 0, 0, 0, 0);
        Particle izq = new Particle(0, p.y, p.z, 0, 0, 0, 0, 0);
        Particle front = new Particle(p.x, p.y, s.dep, 0, 0, 0, 0, 0);
        Particle back = new Particle(p.x, p.y, 0, 0, 0, 0, 0, 0);
        Particle inf = new Particle(p.x, 0, p.z, 0, 0, 0, 0, 0);


        // Comparo contra las paredes ...

        if( p.getOverlap(sup) > 0 ) {
            // Nos traemos FN y FT
            fn = getFN(p, sup, s);
            ft = getFT(p, sup, s);
            ft2 = getFT2(p, sup, s);
            // Calculamos las fuerzas totales en X e Y ...
            double xDiff = sup.x - p.x;
            double yDiff = sup.y - p.y;
            double zDiff = sup.z - p.z;
            double distance = p.getDistance(sup);
            double enx = xDiff/distance;
            double eny = yDiff/distance;
            double enz = zDiff/distance;
            Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
            Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
            Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);
            fxT += fn * enx + ft * (-eny) + ft2 * et2.first;
            fyT += fn * eny + ft * enx + ft2 * et2.second;
            fzT += fn * enz + ft * enz + ft2 * et2.third;
            if(fxT > 1000 || fyT > 1000 || fzT > 1000){
                System.out.println("ME fui a la mierda SUP");
                System.exit(1);
            }
        }
        ////////////////////////////////////////////
        //
        //
        // TODO: REVISAR ESTO FALTA USAR EL FLOORZ
        // TODO: Aca no seria armar por ejemplo el cuadrado y si no esta en cuadrado aplicas fuerza
        // y sino que siga. Y este cuadrado seria usado X e Y, no Z si no entiendo mal.
        // 
        ////////////////////////////////////////////
        // dep
        // ____________________ 
        // |                  |
        // |     ________     | z<dep-fz
        // |    |        |    |
        // |    |        |    |
        // |    |________|    | z>fz
        // |                  |
        // |__________________|
        // 0  x>fx  x<wid-fx   wid
        if(check_outside_square(p, s) && p.getOverlap(inf) > 0 && !p.ignore){

            // if(!(!check_outside_square(s.prevParticles.get(index),s) && p.getOverlap(inf) > 0)){

            // }
            // if(!(p.x >= s.floorX + p.r && p.x <= s.wid - s.floorX - p.r && p.z >= s.floorZ + p.r && p.z <= s.dep - s.floorZ - p.r) && p.getOverlap(inf) > 0 ) {
                // Nos traemos FN y FT
            fn = getFN(p, inf, s);
            ft = getFT(p, inf, s);
            ft2 = getFT2(p, inf, s);
            // Calculamos las fuerzas totales en X e Y ...
            double xDiff = inf.x - p.x;
            double yDiff = inf.y - p.y;
            double zDiff = inf.z - p.z;
            double distance = p.getDistance(inf);
            double enx = xDiff/distance;
            double eny = yDiff/distance;
            double enz = zDiff/distance;
            Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
            Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
            Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);
            fxT += fn * enx + ft * (-eny) + ft2 * et2.first;
            fyT += fn * eny + ft * enx + ft2 * et2.second;
            fzT += fn * enz + ft * enz + ft2 * et2.third;
            if(fxT > 1000 || fyT > 1000 || fzT > 1000){
                System.out.println(fxT + " - " + fyT + " - " + fzT);
                System.out.println(p.x + " - " + p.y + " - " + p.z);
                System.out.println(s.prevParticles.get(index).x + " - " + s.prevParticles.get(index).y + " - " + s.prevParticles.get(index).z);
                System.out.println(p.r);
                System.out.println("ME fui a la mierda INF - " + p.getOverlap(inf));
                System.exit(1);
            }
        } else if(!check_outside_square(p, s) && p.getOverlap(inf) > 0){
            p.ignore = true;
        }

        if( p.getOverlap(izq) > 0 ) {
            // Nos traemos FN y FT
            fn = getFN(p, izq, s);
            ft = getFT(p, izq, s);
            ft2 = getFT2(p, izq, s);
            // Calculamos las fuerzas totales en X e Y ...
            double xDiff = izq.x - p.x;
            double yDiff = izq.y - p.y;
            double zDiff = izq.z - p.z;
            double distance = p.getDistance(izq);
            double enx = xDiff/distance;
            double eny = yDiff/distance;
            double enz = zDiff/distance;
            Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
            Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
            Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);
            fxT += fn * enx + ft * (-eny) + ft2 * et2.first;
            fyT += fn * eny + ft * enx + ft2 * et2.second;
            fzT += fn * enz + ft * enz + ft2 * et2.third;
            if(fxT > 1000 || fyT > 1000 || fzT > 1000){
                System.out.println("ME fui a la mierda IZQ");
                System.exit(1);
            }
        }

        if( p.getOverlap(der) > 0 ) {
            // Nos traemos FN y FT
            fn = getFN(p, der, s);
            ft = getFT(p, der, s);
            ft2 = getFT2(p, der, s);
            // Calculamos las fuerzas totales en X e Y ...
            double xDiff = der.x - p.x;
            double yDiff = der.y - p.y;
            double zDiff = der.z - p.z;
            double distance = p.getDistance(der);
            double enx = xDiff/distance;
            double eny = yDiff/distance;
            double enz = zDiff/distance;
            Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
            Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
            Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);
            fxT += fn * enx + ft * (-eny) + ft2 * et2.first;
            fyT += fn * eny + ft * enx + ft2 * et2.second;
            fzT += fn * enz + ft * enz + ft2 * et2.third;
            if(fxT > 1000 || fyT > 1000 || fzT > 1000){
                System.out.println("ME fui a la mierda DER");
                System.exit(1);
            }
        }
        if( p.getOverlap(front) > 0 ) {
            // Nos traemos FN y FT
            fn = getFN(p, front, s);
            ft = getFT(p, front, s);
            ft2 = getFT2(p, front, s);
            // Calculamos las fuerzas totales en X e Y ...
            double xDiff = front.x - p.x;
            double yDiff = front.y - p.y;
            double zDiff = front.z - p.z;
            double distance = p.getDistance(front);
            double enx = xDiff/distance;
            double eny = yDiff/distance;
            double enz = zDiff/distance;
            Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
            Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
            Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);
            fxT += fn * enx + ft * (-eny) + ft2 * et2.first;
            fyT += fn * eny + ft * enx + ft2 * et2.second;
            fzT += fn * enz + ft * enz + ft2 * et2.third;
            if(fxT > 1000 || fyT > 1000 || fzT > 1000){
                System.out.println("ME fui a la mierda FRONT");
                System.exit(1);
            }
        }
        if( p.getOverlap(back) > 0 ) {
            // Nos traemos FN y FT
            fn = getFN(p, back, s);
            ft = getFT(p, back, s);
            ft2 = getFT2(p, back, s);
            // Calculamos las fuerzas totales en X e Y ...
            double xDiff = back.x - p.x;
            double yDiff = back.y - p.y;
            double zDiff = back.z - p.z;
            double distance = p.getDistance(back);
            double enx = xDiff/distance;
            double eny = yDiff/distance;
            double enz = zDiff/distance;
            Triple<Double, Double, Double> en = new Triple<Double, Double, Double>(enx, eny, enz);
            Triple<Double, Double, Double> et = new Triple<Double, Double, Double>(-eny, enx, enz);
            Triple<Double, Double, Double> et2 = getPerpendicularVersor(en, et);
            fxT += fn * enx + ft * (-eny) + ft2 * et2.first;
            fyT += fn * eny + ft * enx + ft2 * et2.second;
            fzT += fn * enz + ft * enz + ft2 * et2.third;
            if(fxT > 1000 || fyT > 1000 || fzT > 1000){
                System.out.println("ME fui a la mierda BACK");
                System.exit(1);
            }
        }
        fyT -= 9.8 * (p.m);
        return new Triple<Double,Double, Double>(fxT, fyT, fzT);
    }
    
    static boolean check_outside_square(Particle p, Sylo s){
        return (p.x <= s.floorX || p.x >= s.wid - s.floorX) && (p.z <= s.floorZ || p.z  >= s.dep - s.floorZ);
    }
}
