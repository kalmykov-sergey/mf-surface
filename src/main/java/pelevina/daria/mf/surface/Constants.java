package pelevina.daria.mf.surface;

import org.apfloat.Apfloat;

public interface Constants {
    int N = 200;
    int H_0 = 350;
    double R = 0.06647221276;
    double left = R; // 0.2104708820051187;
    double right = 5*R; //0.2104708820070533;
    double alpha = -0.5410520681;
    double Eps = 0.01;
    Apfloat xi_inf = new Apfloat(0.004973591970);
    double P1 = 290.3761507;
    double G = 697.8648613;
    Apfloat Ms = new Apfloat(9.6);
    Apfloat C = new Apfloat(0.134123123123);
    double CONTAINER_RADIUS = 1.0;
    double CONTAINER_HEIGHT = 0.5;
}
