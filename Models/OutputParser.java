package Models;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class OutputParser {

    private static boolean firstEj1 = true;
    private static final String xyz_fn = "XYZ/output.xyz";

    public static void writeUniverse(List<Particle> particles, List<Particle> borders, double t) {
        try {
            StringBuilder dump = new StringBuilder(particles.size() + 8 + "\n" + "Time=" + t + "\n");
            for (Particle p : borders) {
                dump.append(200).append(" ");
                dump.append(p.x).append(" ")
                        .append(p.y).append(" ")
                        .append(p.z).append(" ")
                        .append(p.getVelocity()).append(" ")
                        .append(p.r).append("\n");
            }
            for (Particle p : particles) {
                dump.append(200).append(" ");
                dump.append(p.x).append(" ")
                        .append(p.y).append(" ")
                        .append(p.z).append(" ")
                        .append(p.getVelocity()).append(" ")
                        .append(p.r).append("\n");
            }
            appendToEndOfFile(xyz_fn,dump.toString());
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    public static void parseEj1_2(double dt, int count, String filename){
        try {
            StringBuilder dump = new StringBuilder("");
            if(firstEj1 == true){
                dump.append("Dt,N\n");
                firstEj1 = false;
            }
            dump.append(dt).append(",").append(count).append("\n");
            appendToEndOfFile(filename,dump.toString());
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    public static void parseEj3_4(double dt, double ke, String filename){
        try {
            StringBuilder dump = new StringBuilder("");
            if(firstEj1 == true){
                dump.append("Dt,KE\n");
                firstEj1 = false;
            }
            dump.append(dt).append(",").append(ke).append("\n");
            appendToEndOfFile(filename,dump.toString());
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    private static void appendToEndOfFile(String file,String text) throws IOException {
        FileWriter fw = new FileWriter(file, true);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write(text);
        bw.close();
    }

    public static void createCleanUniverseFile(String fn) {
        Path fileToDeletePath = Paths.get(fn);
        try {
            Files.deleteIfExists(fileToDeletePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void createCleanCSVFile(String fn) {
        Path fileToDeletePath = Paths.get(fn);
        firstEj1 = true;
        try {
            Files.deleteIfExists(fileToDeletePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
