package pelevina.daria.mf.surface;


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
        double rho;
        double u;
        double teta;

        double getX() {
            return rho;
        }

        double getY() {
            return u;
        }

        double getAngle() {
            return teta;
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
