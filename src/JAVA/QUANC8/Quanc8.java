package JAVA.QUANC8;

import java.util.function.Function;

public class Quanc8 {

    private Function<Double, Double> FUN;
    private double A;
    private double B;
    private double ABSERR;
    private double RELERR;
    private double RESULT;
    private double ERREST;
    private int NOFUN;
    private double FLAG;

    private double[] QRIGHT;
    private double[] F;
    private double[] X;
    private double[][] FSAVE;
    private double[][] XSAVE;

    private int LEVMIN, LEVMAX, LEVOUT, NOMAX, NOFIN, LEV, NIM, J, I;
    private double W0;
    private double W1;
    private double W2;
    private double W3;
    private double W4;
    private double COR11;
    private double AREA;
    private double X0;

    private double F0;
    private double STONE;
    private double STEP;
    private double QLEFT, QNOW, QDIFF, QPREV, TOLERR, ESTERR;

    private boolean calculated = false;

    public Quanc8(Function<Double, Double> FUN, double a, double b, double ABSERR, double RELERR) {
        this.FUN = FUN;
        A = a;
        B = b;
        this.ABSERR = ABSERR;
        this.RELERR = RELERR;
    }

    public boolean isCalculated() {
        return calculated;
    }

    public double getRESULT() {
        if(!calculated) {
            throw new RuntimeException("The QUANC8 has not been counted yet, call the function calculate().");
        }
        return RESULT;
    }

    public double getERREST() {
        if(!calculated) {
            throw new RuntimeException("The QUANC8 has not been counted yet, call the function calculate().");
        }
        return ERREST;
    }

    public int getNOFUN() {
        if(!calculated) {
            throw new RuntimeException("The QUANC8 has not been counted yet, call the function calculate().");
        }
        return NOFUN;
    }

    public double getFLAG() {
        if(!calculated) {
            throw new RuntimeException("The QUANC8 has not been counted yet, call the function calculate().");
        }
        return FLAG;
    }

    public void calculate() {
        if(calculated) {
            throw new RuntimeException("Already calculated!");
        }
        QRIGHT = new double[32];
        F = new double[17];
        X = new double[17];
        FSAVE = new double[9][31];
        XSAVE = new double[9][31];

        LEVMIN = 1;
        LEVMAX = 30;
        LEVOUT = 6;
        NOMAX = 5000;
        NOFIN = (int) (NOMAX - (8 * (LEVMAX - LEVOUT + Math.pow(2, (double) (LEVOUT + 1)))));

        W0 = 3956.0 / 14175.0;
        W1 = 23552.0 / 14175.0;
        W2 = -3712.0 / 14175.0;
        W3 = 41984.0 / 14175.0;
        W4 = -18160.0 / 14175.0;

        FLAG = 0.0;
        RESULT = 0.0;
        COR11 = 0.0;
        ERREST = 0.0;
        AREA = 0.0;
        NOFUN = 0;

        if(A == B) {
            throw new IllegalArgumentException("Borders can't be equal");
        }
        LEV = 0;
        NIM = 1;
        X0 = A;
        X[16] = B;
        QPREV = 0.0;
        F0 = FUN.apply(X0);
        STONE = (B - A) / 16.0;
        X[8] = (X0 + X[16]) / 2.0;
        X[4] = (X0 + X[8]) / 2.0;
        X[12] = (X[8] + X[16]) / 2.0;
        X[2] = (X0 + X[4]) / 2.0;
        X[6] = (X[4] + X[8]) / 2.0;
        X[10] = (X[8] + X[12]) / 2.0;
        X[14] = (X[12] + X[16]) / 2.0;
        for (J = 2; J <= 16; J = J + 2) {
            F[J] = FUN.apply(X[J]);
        }
        NOFUN = 9;
        trenta();
    }



    private void trenta() {
        X[1] = (X0 + X[2]) / 2.0;
        F[1] = FUN.apply(X[1]);
        for (J = 3; J <= 15; J = J + 2) {
            X[J] = (X[J - 1] + X[J + 1]) / 2.0;
            F[J] = FUN.apply(X[J]);
        }
        NOFUN = NOFUN + 8;
        STEP = (X[16] - X0) / 16.0;
        QLEFT = (W0 * (F0 + F[8]) + W1 * (F[1] + F[7]) + W2 * (F[2] + F[6]) + W3 * (F[3] + F[5])
                + W4 * F[4]) * STEP;
        QRIGHT[LEV + 1] = (W0 * (F[8] + F[16]) + W1 * (F[9] + F[15]) + W2 * (F[10] + F[14])
                + W3 * (F[11] + F[13]) + W4 * F[12]) * STEP;
        QNOW = QLEFT + QRIGHT[LEV + 1];
        QDIFF = QNOW - QPREV;
        AREA = AREA + QDIFF;
        ESTERR = Math.abs(QDIFF) / 1023.0;
        double tolerr = (RELERR * Math.abs(AREA)) * (STEP / STONE);
        if (ABSERR > tolerr)
            TOLERR = ABSERR;
        else
            TOLERR = tolerr;
        if (LEV < LEVMIN) {
            cinquanta();
        }
        else if (LEV >= LEVMAX) {
            sessantadue();
        }
        else if (NOFUN > NOFIN) {
            sessanta();
        }
        else if (ESTERR <= TOLERR) {
            settanta();
        } else {
            cinquanta();
        }
    }
    private void cinquanta() {
        NIM = 2 * NIM;
        LEV = LEV + 1;


        for (I = 1; I <= 8; I++) {
            FSAVE[I][LEV] = F[I + 8];
            XSAVE[I][LEV] = X[I + 8];
        }
        QPREV = QLEFT;
        for (I = 1; I <= 8; I++) {
            J = -I;
            F[2 * J + 18] = F[J + 9];
            X[2 * J + 18] = X[J + 9];
        }
        trenta();
    }
    private void sessanta() {
        NOFIN = 2 * NOFIN;
        LEVMAX = LEVOUT;
        FLAG = FLAG + ((B - X0) / (B - A));
        settanta();
    }
    private void sessantadue() {
        FLAG = FLAG + 1.0;
        settanta();
    }
    private void settanta() {
        RESULT = RESULT + QNOW;
        ERREST = ERREST + ESTERR;
        COR11 = COR11 + QDIFF / 1023.0;

        while (NIM % 2 != 0) {
            NIM = NIM / 2;
            LEV = LEV - 1;
        }
        NIM = NIM + 1;
        if (LEV <= 0) {
            ottanta();
            return;
        }
        QPREV = QRIGHT[LEV];
        X0 = X[16];
        F0 = F[16];
        for (I = 1; I <= 8; I++) {
            F[2 * I] = FSAVE[I][LEV];
            X[2 * I] = XSAVE[I][LEV];
        }
        trenta();
    }
    private void ottanta() {
        calculated = true;
        RESULT = RESULT + COR11;
        if (ERREST == 0.0)
            return;
        while (Math.abs(RESULT) + (ERREST) == Math.abs(RESULT))
            ERREST = 2.0 * (ERREST);
    }

    @Override
    public String toString() {
        if(!calculated) {
            return "Quanc8{" +
                    "A=" + A +
                    ", B=" + B +
                    ", ABSERR=" + ABSERR +
                    ", RELERR=" + RELERR;
        }
        return "Quanc8{" +
                "A=" + A +
                ", B=" + B +
                ", ABSERR=" + ABSERR +
                ", RELERR=" + RELERR +
                ", RESULT=" + RESULT +
                ", ERREST=" + ERREST +
                ", NOFUN=" + NOFUN +
                ", FLAG=" + FLAG +
                '}';
    }
}
