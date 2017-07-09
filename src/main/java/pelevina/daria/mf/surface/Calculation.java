package pelevina.daria.mf.surface;

import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;

import static java.lang.Math.*;
import static pelevina.daria.mf.surface.Constants.*;
import static pelevina.daria.mf.surface.CrossingFlag.*;

public class Calculation {

    private Apfloat FIVE = new Apfloat(5);
    private Data data;
    private CrossingFlag flag;
    private double initU0;
    private double leftBorder;
    private double rightBorder;

    public Calculation(Data data) {
        this.data = data;
    }

    public void perform() {
        leftBorder = R.doubleValue();
        rightBorder = R.multiply(FIVE).doubleValue();
        flag = NoCrossing;
        flag = oneShot();
    }

    public void next() {
        if (shootFinished()) {
            System.out.println("Here must be volume calc...");
            return;
        }
        flag = oneShot();
    }

    private boolean shootFinished() {
        switch (flag) {
            case WallCrossing:
                double tetac = data.points[data.lastIndex].teta;
                System.out.println("tetac = " + tetac);
                System.out.println("alpha = " + Constants.alpha);
                if (tetac - Constants.alpha.doubleValue()  < -Eps) {
                    rightBorder = initU0;
                } else if (tetac - Constants.alpha.doubleValue()  > Eps) {
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
                if (data.points[data.lastIndex].u < initU0) {
                    leftBorder = initU0;
                } else {
                    rightBorder = initU0;
                }
                break;
        }
        return false;
    }

    private CrossingFlag oneShot() {
        initU0 = (leftBorder + rightBorder) / 2;
        System.out.println("=====================");
        System.out.println("left = " + leftBorder);
        System.out.println("right = " + rightBorder);
        System.out.println("u0 = " + initU0);
        Data.Point initPoint = data.points[0];
        initPoint.u = initU0;
        initPoint.rho = 0;
        initPoint.teta = 0;

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
        double k13 = h * f3(previous.rho, previous.u, Constants.C.doubleValue(), Constants.H_0, previous.teta);
        double k21 = h * f1(previous.teta + k13);
        double k22 = h * f2(previous.teta + k13);
        double k23 = h * f3(previous.rho + k11, previous.u + k12, Constants.C.doubleValue(), Constants.H_0, previous.teta + k13);
        current.rho = previous.rho + (k11 + k21) * 0.5;
        current.u = previous.u + (k12 + k22) * 0.5;
        current.teta = previous.teta + (k13 + k23) * 0.5;
    }

    private CrossingFlag checkCrossing(CrossingFlag previousFlag, int currentPointIndex) {
        Data.Point current = data.points[currentPointIndex];
        CrossingFlag flag = NoCrossing;
        // CHECKING CROSS //
        if ((current.rho < 0) && (current.u < initU0))
            flag = ClockWise;
        if ((current.rho < 0) && (current.u > initU0))
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
        if (current.u < -Constants.R.doubleValue()) {
            flag = BottomCrossing;
        }
        if (current.rho > Constants.CONTAINER_RADIUS) {
            flag = WallCrossing;
        }

        if (pow(current.rho, 2) + pow(current.u, 2) < ApfloatMath.pow(R, 2).doubleValue()) {
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
            double f1 = (current.rho - prevLast.rho) / (last.rho - prevLast.rho) - (current.u - prevLast.u) / (last.u - prevLast.u);
            double f2 = (prevCurrent.rho - prevLast.rho) / (last.rho - prevLast.rho) - (prevCurrent.u - prevLast.u) / (last.u - prevLast.u);
            double g1 = (last.rho - prevCurrent.rho) / (current.rho - prevCurrent.rho) - (last.u - prevCurrent.u) / (current.u - prevCurrent.u);
            double g2 = (prevLast.rho - prevCurrent.rho) / (current.rho - prevCurrent.rho) - (prevLast.u - prevCurrent.u) / (current.u - prevCurrent.u);
            double k = (current.rho - last.rho) * (prevLast.u - current.u) - (current.u - last.u) * (prevLast.rho - current.rho);
            if ((f1 * f2 < 0) && (g2 * g1 < 0)) {
                if (k > 0) return CounterClockWise;
                if (k < 0) return ClockWise;
            }
        }
        return NoCrossing;
    }

    private static double H(double r, double z) {
        double A = - ApfloatMath.pow(R, 3).doubleValue();
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

    private static Apfloat f1(Apfloat x) {
        return ApfloatMath.cos(x);
    }

    private static double f2(double y) {
        return -sin(y);
    }

    private static Apfloat f2(Apfloat x) {
        return ApfloatMath.sin(x).negate();
    }

}
