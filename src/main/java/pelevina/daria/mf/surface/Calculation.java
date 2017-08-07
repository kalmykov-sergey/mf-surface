package pelevina.daria.mf.surface;

import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.apfloat.Apint;

import static java.lang.Math.*;
import static pelevina.daria.mf.surface.Constants.*;
import static pelevina.daria.mf.surface.CrossingFlag.*;

public class Calculation {
    private Data data;
    private CrossingFlag flag;
    private Apfloat initU0;
    private Apfloat leftBorder;
    private Apfloat rightBorder;

    public Calculation(Data data) {
        this.data = data;
        leftBorder = new Apfloat(Constants.left);
        rightBorder = new Apfloat(Constants.right);
        flag = NoCrossing;
    }

    public void perform() {

        flag = oneShot();
    }

    public void next() {
        flag = oneShot();
        if (shootFinished()) {
            System.out.println("Here must be volume calc.??..");
            volumeCalculate();
            System.exit(0);
        }
    }

    public boolean pause() {
        boolean pause = !CrossingFlag.CONTINUE.contains(flag);
        if (pause) {
            System.out.println("\n ==== PAUSE ==== \n");
        }
        return false;
        //return pause;
    }

    private boolean shootFinished() {
        switch (flag) {
            case WallCrossing:
                double tetac = data.points[data.lastIndex].teta.doubleValue();
                System.out.println("tetac = " + tetac);
                System.out.println("alpha = " + Constants.alpha);
                if (tetac - Constants.alpha < -Eps) {
                    rightBorder = initU0;
                } else if (tetac - Constants.alpha > Eps) {
                    leftBorder = initU0;
                } else {
                    return true;
                }
                break;
            case ClockWise:
            case SphereCrossing:
            case BottomCrossing:
                leftBorder = initU0;
                break;
            case CounterClockWise:
                rightBorder = initU0;
                break;
            case GoesToInfty:
            case NoCrossing:
                if (data.points[data.lastIndex].u.compareTo(initU0) < 0) {
                    leftBorder = initU0;
                } else {
                    rightBorder = initU0;
                }
                break;
        }
        return false;
    }

    private void volumeCalculate() {
        double tetac = data.points[data.lastIndex].teta.doubleValue();
        System.out.println("Quitting: TETAsm=" + tetac);
        System.out.println(String.format("h0=%f, C=%f\n", initU0.doubleValue(), C.doubleValue()));

        // CALCULATIG OF THE VOLUME //

        double Vol = getVolume();

        Data.Point touch = getTouchPoint();
        // System.out.println(String.format("free_mag_energy=%f; gravity energy=%f; surftension energy=%f\n",fme, fge, fse);
        System.out.println(String.format("V0l=%f, rho0=%f\n", Vol, touch.rho.doubleValue()));
        for (int i = 0; i < data.lastIndex; i++) {
            Data.Point current = data.points[i];
            System.out.println("["
                    + String.valueOf(current.rho.doubleValue()).replace(',', '.')
                    + ", "
                    + String.valueOf(current.u.doubleValue()).replace(',', '.')
                    + "], ");
        }
    }

    private Data.Point getTouchPoint() {
        Data.Point touch = data.points[0];
        double delta = 1;
        for (int i = 0; i < data.lastIndex + 1; i++) {
            Data.Point current = data.points[i];
            if (current.R2() - R * R < delta) {
                delta = current.R2() - R * R;
                touch = current;
            }
        }
        return touch;
    }

    private double getVolume() {
        double vol = 0;
        for (int i = 1; i < data.lastIndex; i++) {
            Data.Point current = data.points[i];
            Data.Point previous = data.points[i - 1];
            vol = vol + 0.5 * Math.PI * (previous.u.doubleValue()  + current.u.doubleValue() ) * (pow(current.rho.doubleValue(), 2) - pow(previous.rho.doubleValue(), 2));
        }
        //vol = vol - 4 * Math.PI * R * R * R / 3;
        return vol;
    }

    private CrossingFlag oneShot() {
        initU0 = leftBorder.add(rightBorder).divide(new Apfloat(2));
        System.out.println("==== calculation started ======");
        System.out.println("left = " + leftBorder);
        System.out.println("right = " + rightBorder);
        if (leftBorder.equals(rightBorder)) {
            System.out.println("too few digits");
            System.exit(1);
        }
        System.out.println("u0 = " + initU0);
        Data.Point initPoint = data.points[0];
        initPoint.u = initU0;
        initPoint.rho = Apfloat.ZERO;
        initPoint.teta = Apfloat.ZERO;

        for (int currentPointIndex = 1; currentPointIndex < data.points.length; currentPointIndex++) {
            Data.Point current = data.points[currentPointIndex];
            calcStep(data.points[currentPointIndex], data.points[currentPointIndex - 1]);
            flag = checkCrossing(currentPointIndex);
            if (flag != NoCrossing) {
                data.lastIndex = currentPointIndex;
                System.out.println("===== calculation finished ======");
                return flag;
            }

        }
        System.out.println("===== calculation finished ======");
        System.out.println("end of points");
        return GoesToInfty;
    }

    private void calcStep(Data.Point current, Data.Point previous) {
        Apfloat h = new Apfloat(1.0).divide(new Apint(N));
        Apfloat k11 = h.multiply(f1(previous.teta));
        Apfloat k12 = h.multiply(f2(previous.teta));
        Apfloat k13 = h.multiply(f3(previous.rho, previous.u, Constants.C.doubleValue(), Constants.H_0, previous.teta));
        Apfloat k21 = h.multiply(f1(previous.teta.add(k13)));
        Apfloat k22 = h.multiply(f2(previous.teta.add(k13)));
        Apfloat k23 = h.multiply(f3(previous.rho.add(k11), previous.u.add(k12), Constants.C.doubleValue(), Constants.H_0, previous.teta.add(k13)));
        current.rho = previous.rho.add((k11.add(k21)).divide(new Apfloat(2)));
        current.u = previous.u.add((k12.add(k22)).divide(new Apfloat(2)));
        current.teta = previous.teta.add((k13.add(k23)).divide(new Apint(2)));
    }

    private CrossingFlag checkCrossing(int currentPointIndex) {
        Data.Point current = data.points[currentPointIndex];
        CrossingFlag flag = NoCrossing;
        // CHECKING CROSS //
        if ((current.rho.compareTo(Apfloat.ZERO) < 0) && (current.u.compareTo(initU0) < 0))
            flag = ClockWise;
        if ((current.rho.compareTo(Apfloat.ZERO) < 0) && (current.u.compareTo(initU0) > 0))
            flag = CounterClockWise;
        if (currentPointIndex > 2) {
            CrossingFlag crtmp = cross(currentPointIndex);
            if (
                    ((crtmp == CounterClockWise) && (flag == ClockWise))
                            ||
                            ((crtmp == ClockWise) && (flag == CounterClockWise))
                    )
                flag = NoCrossing;
            else if (crtmp != NoCrossing)
                flag = crtmp;
        }
        if (current.u.compareTo(new Apfloat(-Constants.R)) < 0) {
            flag = BottomCrossing;
        }
        if (current.rho.compareTo(new Apfloat(Constants.CONTAINER_RADIUS)) > 0) {
            flag = WallCrossing;
        }

        if (pow(current.rho.doubleValue(), 2) + pow(current.u.doubleValue(), 2) < pow(R, 2)) {
            System.out.println(current);
            flag = SphereCrossing;
        }

        return flag;
    }

    // CHECKING: CURVE CROSS HIMSELF AND DIRECTION OF THE TURN //
    private CrossingFlag cross(int currentPointIndex) {
        for (int j = 2; j < currentPointIndex - 1; j++) {
            Data.Point last = data.points[currentPointIndex];
            Data.Point prevLast = data.points[currentPointIndex - 1];
            Data.Point current = data.points[j];
            Data.Point prevCurrent = data.points[j - 1];
            Apfloat f1 = crossingDiff(current, last, prevLast);
            Apfloat f2 = crossingDiff(prevCurrent, last, prevLast);
            Apfloat g1 = crossingDiff(last, current, prevCurrent);
            Apfloat g2 = crossingDiff(prevLast, current, prevCurrent);
            if ((f1.multiply(f2).compareTo(Apfloat.ZERO) < 0) && (g2.multiply(g1).compareTo(Apfloat.ZERO) < 0)) {
                double k = (current.rho.subtract(last.rho).doubleValue()) * (prevLast.u.subtract(current.u).doubleValue()) - (current.u.subtract(last.u).doubleValue()) * (prevLast.rho.subtract(current.rho).doubleValue());
                if (k > 0) return CounterClockWise;
                if (k < 0) return ClockWise;
            }
        }
        return NoCrossing;
    }

    private Apfloat crossingDiff(Data.Point first, Data.Point second, Data.Point base) {
        return (
                (first.rho.subtract(base.rho).divide(second.rho.subtract(base.rho)))
                        .subtract((first.u.subtract(base.u).divide(second.u.subtract(base.u))))
        );
    }

    private static double H(double r, double z) {
        double A = -pow(R, 3);
        double z2 = z * z;
        double r2 = r * r;
        double sum = z2 + r2;
        double H2 = 9 * A * A * z2 * r2 / pow(sum, 5)
                + (1 + A * (1 / pow(sum, 1.5) - 3 * z2 / pow(sum, 2.5)))
                * (1 + A * (1 / pow(sum, 1.5) - 3 * z2 / pow(sum, 2.5)));
        double H = sqrt(H2);
//        System.out.println("H=" + H);
        return H;
    }

    private static double langeven(double r, double z, double H0) {
        double xi = Constants.xi_inf.doubleValue() * H0;
        double H_nd = H(r, z);
        double P = ((exp(xi * H_nd) + exp(-xi * H_nd))
                / (exp(xi * H_nd) - exp(-xi * H_nd)))
                - 1 / (xi * H_nd);
        System.out.println("P=" + P);
        return Constants.Ms.doubleValue() * P;
    }

    private static Apfloat P(double r, double z, double H0) {
        double xi = Constants.xi_inf.doubleValue() * H0;
        double H_nd = H(r, z);
        double pp = log(exp(xi * H_nd) / 2 - exp(-xi * H_nd) / 2) - log(xi * H_nd);
        //printf("P=%f\n",pp);
        return new Apfloat(pp);
    }

    private static Apfloat f3(Apfloat r, Apfloat u, double C, double H0, Apfloat theta) {
        Apfloat k = (u.multiply(new Apfloat(Constants.G)).negate()).add(new Apfloat(Constants.P1).multiply(P(r.doubleValue(), u.doubleValue(), H0)));
        //printf("-Ah+P=%f,C-Ah+P=%f  ",k,k+C);
        if (abs(r.doubleValue()) < Eps) {
            return ((k.add(new Apfloat(C))).divide(new Apint(2)));
        }
        return k.add(new Apfloat(C)).subtract(ApfloatMath.sin(theta).divide(r));
    }

    private static Apfloat f1(Apfloat x) {
        return ApfloatMath.cos(x);
    }

    private static Apfloat f2(Apfloat x) {
        return ApfloatMath.sin(x).negate();
    }

}
