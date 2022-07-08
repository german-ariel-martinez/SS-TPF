import Models.Triple;
import Models.*;

public class App {
    public static void main(String[] args) {
        double dt = 0.00001;
        Sylo sylo = new Sylo(1, 0.3, 0.3, 0.17, dt, 2);
        OutputParser.createCleanUniverseFile("XYZ/output.xyz");
        System.out.println("Entro al populate\n");
        sylo.populate(0.005);
        sylo.simulateUniverse(50);
    }
}