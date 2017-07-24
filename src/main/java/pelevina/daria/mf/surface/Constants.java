package pelevina.daria.mf.surface;

import org.apfloat.Apfloat;

public interface Constants {
    int N = 200;
    int H_0 = 1;
    double R = 0.9756097561;
    double left = R; // 0.2104708820051187;
    double right = 5*R; //0.2104708820070533;
    double alpha = 0.35;
    double Eps = 0.01;
    Apfloat xi_inf = new Apfloat(0.007587829123);
    double P1 = 48.63051;
    double G = 3.465675675;
    Apfloat Ms = new Apfloat(12);
    Apfloat C = new Apfloat(20.38632750);
    double CONTAINER_RADIUS = 1.0;
    double CONTAINER_HEIGHT = 1.5;
}
