package pelevina.daria.mf.surface;


import org.apfloat.Apfloat;

public class Data {

    final Point[] points;
    public int lastIndex;

    public Data(int numberOfPoints) {
        this.points = new Point[numberOfPoints];
        for (int i = 0; i < points.length; i++) {
            points[i] = new Point();
        }
        lastIndex = this.points.length - 1;
    }

    public static class Point {
        Apfloat rho;
        Apfloat u;
        Apfloat teta;

        double getX() {
            return rho.doubleValue();
        }

        double getY() {
            return u.doubleValue();
        }

        double getAngle() {
            return teta.doubleValue();
        }

        @Override
        public String toString() {
            return "Point{" +
                    "rho=" + getX() +
                    ", u=" + getY() +
                    ", teta=" + getAngle() +
                    '}';
        }
    }
}
