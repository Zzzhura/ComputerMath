package JAVA.RKF45;


import static java.lang.Math.*;
import static java.lang.System.exit;

public class RKF45 {

    private static final double REMIN = 1.0e-12;
    private static final double EPSILON = 2.2e-16;

    private static int HFAILD, OUTPUT;
    private static double A, AE, DT, EE, EEOET, ESTTOL, ET, HMIN, RER, S,
            SCALE, TOL, TOLN, U26, YPK;
    private static int K, MFLAG;

    private static double[] F1, F2, F3, F4, F5;
    private static double SAVRE, SAVAE;
    private static int KOP, INIT, JFLAG, KFLAG;

    private FuncInterf funcInterf;
    private int NEQN;
    private double[] Y;
    private double[] YP;
    private double T;
    private double TOUT;
    private double RELERR;
    private double ABSERR;
    private double H;
    private int NFE;
    private int MAXNFE;
    private int IFLAG;

    public double getT() {
        return T;
    }

    public double getRELERR() {
        return RELERR;
    }

    public double getH() {
        return H;
    }

    public int getNFE() {
        return NFE;
    }

    public int getIFLAG() {
        return IFLAG;
    }

    public RKF45(FuncInterf funcInterf, int NEQN, double[] y, double[] YP, double t, double TOUT, double RELERR,
                 double ABSERR, double h, int NFE, int MAXNFE, int IFLAG) {
        this.funcInterf = funcInterf;
        this.NEQN = NEQN;
        Y = y;
        this.YP = YP;
        T = t;
        this.TOUT = TOUT;
        this.RELERR = RELERR;
        this.ABSERR = ABSERR;
        H = h;
        this.NFE = NFE;
        this.MAXNFE = MAXNFE;
        this.IFLAG = IFLAG;

        F1 = new double[NEQN];
        F2 = new double[NEQN];
        F3 = new double[NEQN];
        F4 = new double[NEQN];
        F5 = new double[NEQN];
    }

    private void init() throws Exception {
        if (NEQN < 1) {
            IFLAG = 8;
            throw new Exception();
        }
        if ((RELERR < 0.0) || (ABSERR < 0.0)) {
            IFLAG = 8;
            throw new Exception();
        }
        MFLAG = abs(IFLAG);
        if ((MFLAG == 0) || (MFLAG > 8)) {
            IFLAG = 8;
            throw new Exception();
        }

        if (MFLAG != 1) {
            /* NOT the first call: Check continuation possibilities */
            if ((T == TOUT) && (KFLAG != 3)) {
                IFLAG = 8;
                throw new Exception();
            }
            if (MFLAG == 2) {  /*  IFLAG = +2 or -2 */
                switch (KFLAG) {
                    case 3:
                        if (INIT == 0) { /* Reset flag value from previous call */
                            IFLAG = JFLAG;
                            MFLAG = abs(IFLAG);
                        }
                        break;
                    case 4:
                        NFE = 0;  /* Reset function evaluation counter */
                        break;
                    case 5:
                        if (ABSERR == 0.0) exit(0); /* stop here, user did */
                        break;                      /* not heed warning    */
                    case 6:
                        if ((RELERR <= SAVRE) && (ABSERR <= SAVAE)) exit(0);
                        break;                      /* as for case 5       */
                }  /* end switch (*KFLAG) ... */
            } else {  /* IFLAG = 3,4,5,6,7 OR 8 */
                switch (IFLAG) {
                    case 3:
                        IFLAG = JFLAG;
                        if (KFLAG == 3) MFLAG = abs(IFLAG);
                        break;
                    case 4:
                        NFE = 0;
                        IFLAG = JFLAG;
                        if (KFLAG == 3) MFLAG = abs(IFLAG);
                        break;
                    case 5:
                        if (ABSERR > 0.0) {  /* Reset flag from previous call */
                            IFLAG = JFLAG;
                            if (KFLAG == 3) MFLAG = abs(IFLAG);
                        }
                        break;
                    default:
                        exit(0);    /* stop here as the user did not
                                     fix the problem pertaining to
                                     IFLAG = 5, 6, 7 or 8   */
                }  /* end switch (*IFLAG) */
            }  /* end if (MFLAG == 2) ... */
        } /* end if (MFLAG != 1) ... */


 /* -----------------
    Do some real work ...
    -----------------     */

    /* Save input IFLAG and set continuation flag value for subsequent
       input checking */
        JFLAG = IFLAG;
        KFLAG = 0;

        /* Save RELERR and ABSERR for checking input on subsequent calls */
        SAVRE = RELERR;
        SAVAE = ABSERR;

    /* Restrict relative error tolerance to be at least as large as
       2*EPSILON+REMIN to avoid limiting precision difficulties arising
       from impossible accuracy requests  */

        RER = 2.0 * EPSILON + REMIN;
        if (RELERR < RER) {  /* Relative error tolerance too small  */
            RELERR = RER;
            IFLAG = 3;
            KFLAG = 3;
            throw new Exception();
        }

        DT = TOUT - T;


    /* --------------
       Initialization ...
       -------------- */

        if (MFLAG == 1) {  /* set initialization completion indicator,INIT
          set indicator for too many output points,KOP
          evaluate initial derivatives
          set counter for function evaluations,NFE  */
            INIT = 0;
            KOP = 0;
            A = T;
            funcInterf.f(NEQN, A, Y, YP);
            NFE = 1;
            if (T == TOUT) {
                IFLAG = 2;
                throw new Exception();
            }
        }

        if ((MFLAG == 1) || (INIT == 0)) {  /* estimate starting stepsize */
            INIT = 1;
            H = abs(DT);
            TOLN = 0.0;
            for (K = 0; K < NEQN; ++K) {
                TOL = RELERR * abs(Y[K]) + ABSERR;
                if (TOL <= 0.0) continue;
                TOLN = TOL;
                YPK = abs(YP[K]);
                /* *H and TOL/YPK always +ve */
                if ((YPK * pow(H, 5.0)) > TOL) {
                    H = pow((TOL / YPK), 0.2);
                }
            }
            if (TOLN <= 0.0) {
                H = 0.0;
            }
            H = max(H, U26 * max(abs(T), abs(DT)));
            JFLAG = (IFLAG > 0) ? 2 : -2;
        }

        /* Set stepsize for integration in the direction from T to TOUT  */
        H = (DT > 0.0) ? abs(H) : -abs(H);

    /* Test to see if RKF45 is being severely impacted by too many
       output points */
        if (abs(H) >= 2.0 * abs(DT)) ++(KOP);
        if (KOP == 100) {  /* Unnecessary frequency of output  */
            KOP = 0;
            IFLAG = 7;
            throw new Exception();
        }

        /* If too close to output point, extrapolate (linearly) and return  */
        if (abs(DT) <= U26 * abs(T)) {
            for (K = 0; K < NEQN; ++K) {
                Y[K] += DT * YP[K];
            }
            A = TOUT;
            funcInterf.f(NEQN, A, Y, YP);
            ++(NFE);
            T = TOUT;
            IFLAG = 2;
            throw new Exception();
        }

        /* Initialize output point indicator  */
        OUTPUT = 0;  /* FALSE */

    /* To avoid premature underflow in the error tolerance function.
       scale the error tolerances  */
        SCALE = 2.0 / RELERR;
        AE = SCALE * ABSERR;


    /* ------------------------
       Step by step integration ...
       ------------------------ */
        takeSomeSteps();
    }

    public void takeSomeSteps() throws Exception {
        HFAILD = 0;   /* .FALSE. */

        /* Set smallest allowable stepsize */
        HMIN = U26 * abs(T);

    /* Adjust stepsize if necessary to hit the output point.
       look ahead two steps to avoid drastic changes in the stepsize and
       thus lessen the impact of output points on the code.  */
        DT = TOUT - T;
        if (abs(DT) < 2.0 * abs(H)) {
            if (abs(DT) <= abs(H)) {  /* The next successful step will complete the integration to the
           output point  */
                OUTPUT = 1;   /* .TRUE. */
                H = DT;
            } else H = 0.5 * DT;
        }


    /* ----------------------------------------
       Core integrator for taking a single step ...
       ----------------------------------------
       The tolerances have been scaled to avoid premature underflow in
       computing the error tolerance function ET.
       To avoid problems with zero crossings,relative error is measured
       using the average of the magnitudes of the solution at the
       beginning and end of a step.
       The error estimate formula has been grouped to control loss of
       significance.
       To distinguish the various arguments, *H is not permitted
       to become smaller than 26 units of roundoff in T.
       Practical limits on the change in the stepsize are enforced to
       smooth the stepsize selection process and to avoid excessive
       chattering on problems having discontinuities.
       To prevent unnecessary failures, the code uses 9/10 the stepsize
       it estimates will succeed.
       After a step failure, the stepsize is not allowed to increase for
       the next attempted step. This makes the code more efficient on
       problems having discontinuities and more effective in general
       since local extrapolation is being used and extra caution seems
       warranted.
    */

    /* Test number of derivative function evaluations.
       If okay,try to advance the integration from t to t+h  */

        takeAStep();
    }

    private void takeAStep() throws Exception {
        if (NFE > MAXNFE) {  /* Too much work */
            IFLAG = 4;
            KFLAG = 4;
            throw new Exception();
        }

        /* Advance an approximate solution over one step of length H  */
        fehl45(T, H, Y, YP, F1, F2, F3, F4, F5, NEQN);
        NFE += 5;

    /* Compute and test allowable tolerances versus local error estimates
       and remove scaling of tolerances. Note that relative error is
       measured with respect to the average of the magnitudes of the
       solution at the beginning and end of the step.  */
        EEOET = 0.0;
        for (K = 0; K < NEQN; ++K) {
            ET = abs(Y[K]) + abs(F1[K]) + AE;
            if (ET <= 0.0) {   /* Inappropriate error tolerance */
                IFLAG = 5;
                throw new Exception();
            }
            EE = abs((-2090.0 * YP[K] + (21970.0 * F3[K] - 15048.0 * F4[K])) +
                    (22528.0 * F2[K] - 27360.0 * F5[K]));
            EEOET = max(EEOET, EE / ET);
        }

        ESTTOL = abs(H) * EEOET * SCALE / 752400.0;

        if (ESTTOL > 1.0) {  /* --- Unsuccessful step ---
       reduce the stepsize , try again
       the decrease is limited to a factor of 1/10  */
            HFAILD = 1;  /* .TRUE.  */
            OUTPUT = 0;  /* .FALSE. */
            S = 0.1;
            /* ESTTOL always +ve */
            if (ESTTOL < 59049.0) S = 0.9 / pow(ESTTOL, 0.2);
            H *= S;
            if (abs(H) <= HMIN) { /* Requested error unattainable at smallest
             allowable stepsize  */
                IFLAG = 6;
                KFLAG = 6;
                throw new Exception();
            }
            takeAStep();
        }

    /* --- successful step ---
       store solution at T+H
       and evaluate derivatives there  */
        T += H;
        for (K = 0; K < NEQN; ++K) {
            Y[K] = F1[K];
        }
        A = T;
        funcInterf.f(NEQN, A, Y, YP);
        ++NFE;


    /* --- Choose next stepsize. ---
       The increase is limited to a factor of 5.
       If step failure has just occurred, next
       stepsize is not allowed to increase.  */
        S = 5.0;
        if (ESTTOL > 1.889568E-4) S = 0.9 / pow(ESTTOL, 0.2);
        if (HFAILD != 0) S = min(S, 1.0);
        H = (H > 0.0) ? abs(max(S * abs(H), HMIN)) : -abs(max(S * abs(H), HMIN));

    /* ----------------------
       End of core integrator ...
       ---------------------- */


        /* --- Should we take another step? --- */

        if (OUTPUT != 0) {  /* Integration successfully completed over the interval */
            T = TOUT;
            IFLAG = 2;
            return;
        }

        if (IFLAG > 0) {
            takeSomeSteps();
        } else {  /* A single step was taken successfully */
            IFLAG = (-2);
            throw new Exception();
        }

    }

    public void complete() throws Exception {
        init();
    }

    public int fehl45(double T, double H, double[] Y, double[] YP, double[] F1, double[] F2, double[] F3, double[] F4,
                      double[] F5, int NEQN) {
        double CH;
        int K;

        CH = H / 4.0;
        for (K = 0; K < NEQN; ++K) {
            F5[K] = Y[K] + CH * YP[K];
        }
        funcInterf.f(NEQN, T + CH, F5, F1);

        CH = 3.0 * H / 32.0;
        for (K = 0; K < NEQN; ++K) {
            F5[K] = Y[K] + CH * (YP[K] + 3.0 * F1[K]);
        }
        funcInterf.f(NEQN, T + 3.0 * H / 8.0, F5, F2);

        CH = H / 2197.0;
        for (K = 0; K < NEQN; ++K) {
            F5[K] = Y[K] + CH * (1932.0 * YP[K] + (7296.0 * F2[K] -
                    7200.0 * F1[K]));
        }
        funcInterf.f(NEQN, T + 12.0 * H / 13.0, F5, F3);

        CH = H / 4104.0;
        for (K = 0; K < NEQN; ++K) {
            F5[K] = Y[K] + CH * ((8341.0 * YP[K] - 845.0 * F3[K]) +
                    (29440.0 * F2[K] - 32832.0 * F1[K]));
        }
        funcInterf.f(NEQN, T + H, F5, F4);
        CH = H / 20520.0;
        for (K = 0; K < NEQN; ++K) {
            F1[K] = Y[K] + CH * ((-6080.0 * YP[K] + (9295.0 * F3[K] -
                    5643.0 * F4[K])) + (41040.0 * F1[K] - 28352.0 * F2[K]));
        }
        funcInterf.f(NEQN, T + H / 2.0, F1, F5);
        /* --- Compute approximate solution at T+H. --- */

        CH = H / 7618050.0;
        for (K = 0; K < NEQN; ++K) {
            F1[K] = Y[K] + CH * ((902880.0 * YP[K] + (3855735.0 * F3[K] -
                    1371249.0 * F4[K])) + (3953664.0 * F2[K] +
                    277020.0 * F5[K]));
        }

        return 0;
    }


}
