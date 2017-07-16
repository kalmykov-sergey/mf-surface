package pelevina.daria.mf.surface;

import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;

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
    }

    public CrossingFlag perform() {
        leftBorder = new Apfloat(R);
        rightBorder = new Apfloat(R * 5);
        flag = NoCrossing;
        flag = oneShot();
        return flag;
    }

    public CrossingFlag next() {
        if (shootFinished()) {
            System.out.println("Here must be volume calc...");
            return null;
        }
        flag = oneShot();
        return flag;
    }

    private boolean shootFinished() {
        switch (flag) {
            case WallCrossing:
                double tetac = data.points[data.lastIndex].teta.doubleValue();
                System.out.println("tetac = " + tetac);
                System.out.println("alpha = " + Constants.alpha);
                if (tetac - Constants.alpha  < -Eps) {
                    rightBorder = initU0;
                } else if (tetac - Constants.alpha  > Eps) {
                    leftBorder = initU0;
                } else {
                    return true;
                }
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

    private CrossingFlag oneShot() {
        initU0 = leftBorder.add(rightBorder).divide(new Apfloat(2));
        System.out.println("=====================");
        System.out.println("left = " + leftBorder);
        System.out.println("right = " + rightBorder);
        System.out.println("u0 = " + initU0);
        Data.Point initPoint = data.points[0];
        initPoint.u = initU0;
        initPoint.rho = Apfloat.ZERO;
        initPoint.teta = Apfloat.ZERO;

        for (int currentPointIndex = 1; currentPointIndex < data.points.length; currentPointIndex++) {
            Data.Point current = data.points[currentPointIndex];
            calcStep(data.points[currentPointIndex], data.points[currentPointIndex - 1]);
            flag = checkCrossing(flag, currentPointIndex);
            if (flag != NoCrossing) {
                data.lastIndex = currentPointIndex;
                return flag;
            }

        }
        System.out.println("end of points");
        return GoesToInfty;
    }

    private void calcStep(Data.Point current, Data.Point previous) {
        double h = 1.0 / N;
        double k11 = h * f1(previous.teta);
        double k12 = h * f2(previous.teta);
        double k13 = h * f3(previous.rho.doubleValue(), previous.u.doubleValue(), Constants.C.doubleValue(), Constants.H_0, previous.teta.doubleValue());
        double k21 = h * f1(previous.teta.doubleValue() + k13);
        double k22 = h * f2(previous.teta.doubleValue() + k13);
        double k23 = h * f3(previous.rho.doubleValue() + k11, previous.u.doubleValue() + k12, Constants.C.doubleValue(), Constants.H_0, previous.teta.doubleValue() + k13);
        current.rho = previous.rho.add(new Apfloat((k11 + k21) * 0.5));
        current.u = previous.u.add(new Apfloat((k12 + k22) * 0.5));
        current.teta = new Apfloat(previous.teta.doubleValue() + (k13 + k23) * 0.5);
    }

    private CrossingFlag checkCrossing(CrossingFlag previousFlag, int currentPointIndex) {
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
        double A = - pow(R, 3);
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

    private static double P(double r, double z, double H0) {
        double xi = Constants.xi_inf.doubleValue() * H0;
        double H_nd = H(r, z);
        double pp = log(exp(xi * H_nd) / 2 - exp(-xi * H_nd) / 2) - log(xi * H_nd);
        //printf("P=%f\n",pp);
        return pp;
    }

    private static double f3(double r, double u, double C, double H0, double theta) {
        double k = -Constants.G.doubleValue() * u + Constants.P1.doubleValue() * P(r, u, H0);
        double ff;
        //printf("-Ah+P=%f,C-Ah+P=%f  ",k,k+C);
        if (r == 0) {
            ff = (C + k) * 0.5;
            return ff;
        }
        ff = (C + k) - 1.0 * sin(theta) / r;
        return ff;
    }

    private static double f1(double y) {
        return cos(y);
    }

    private static double f1(Apfloat x) {
        return ApfloatMath.cos(x).doubleValue();
    }

    private static double f2(double y) {
        return -sin(y);
    }

    private static double f2(Apfloat x) {
        return ApfloatMath.sin(x).negate().doubleValue();
    }

}
