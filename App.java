import Models.Triple;
import Models.*;

public class App {
    public static void main(String[] args) {
        double dt = 0.00001;
        Sylo sylo = new Sylo(1, 0.4, 0.4,0.1, dt, 0.25);
        OutputParser.createCleanUniverseFile("XYZ/output.xyz");
        System.out.println("Entro al populate\n");
        sylo.populate(0.005);
        sylo.simulateUniverse(50);
    }
}