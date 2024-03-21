package JAVA.ZEROIN;

import java.util.function.Function;

import static java.lang.Math.abs;

public class Zeroin {

    private static final double EPSILON = 2.2e-16;

    private Function<Double, Double> F;
    private double left, right;
    private double tol;
    private int flag;

    public Zeroin(Function<Double, Double> f, double left, double right, double tol) {
        F = f;
        this.left = left;
        this.right = right;
        this.tol = tol;
        this.flag = 0;
    }

    public int getFlag() {
        return flag;
    }

    public double complete() {

        double zero, half, one, two, three;
        double a, b, c = 0.0, d = 0.0, e = 0.0;
        double fa, fb, fc = 0.0, tol1;
        double xm, p, q, r, s;
        int zflag, skip, exit, bracket;
        int nstep, i, nseg;
        double factor, x1, x2, dx, f1, f2;

        zero = 0.0;
        half = 0.5;
        one = 1.0;
        two = 2.0;
        three = 3.0;

        /* Initialization */
        zflag = 0;
        exit = 0;
        a = left;
        b = right;
        fa = F.apply(a);
        fb = F.apply(b);

        /* Check constraints */
        if (tol <= zero || left == right) {
            exit = 1;
            zflag = 2;
        }

        if (fa * (fb / abs(fb)) > zero) {
            /* try to bracket a zero ... */
            bracket = 0;

   /* first check the possibility of an even number of zeros
      within the user supplied range */
            nseg = 10;
            dx = (b - a) / nseg;
            x1 = a;
            f1 = fa;
            for (i = 0; i < nseg; ++i) {
                x2 = x1 + dx;
                f2 = F.apply(x2);
                if (f1 * (f2 / abs(f2)) < zero) {
                    /* this segment brackets a zero */
                    bracket = 1;
                    a = x1;
                    fa = f1;
                    b = x2;
                    fb = f2;
                    break;
                }
                x1 = x2;
                f1 = f2;
            }

            if (bracket == 0) {
                /* now try extending the user supplied range ... */
                factor = 1.6;   /* increase the range by this factor */
                nstep = 20;    /* maximum number of steps */
                x1 = a;
                f1 = fa;
                x2 = b;
                f2 = fb;
                for (i = 0; i < nstep; ++i) {
                    /* extend the range in the downhill direction */
                    if (abs(f1) < abs(f2)) {
                        x1 -= (x2 - x1) * factor;
                        f1 = F.apply(x1);
                    } else {
                        x2 += (x2 - x1) * factor;
                        f2 = F.apply(x2);
                    }
                    if (f1 * (f2 / abs(f2)) <= zero) {
                        /* we have bracketed a zero (or odd number of) */
                        bracket = 1;
                        a = x1;
                        fa = f1;
                        b = x2;
                        fb = f2;
                        break;
                    }
                }
            }
            if (bracket == 0) {
      /* we have been unsuccessful in trying to bracket an
         odd number of zeros */
                exit = 1;
                zflag = 1;
            }
        }


        /* Begin step */
        skip = 0;
        while (exit == 0) {

            if (skip == 0) {
                c = a;        /* ensure that the zero is between b and c */
                fc = fa;
                d = b - a;
                e = d;
            }

            if (abs(fc) < abs(fb)) {
                a = b;             /* swap b and c to give fc >= fb */
                b = c;             /* b is then the best estimate for the zero */
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            /* Convergence test */
            tol1 = two * EPSILON * abs(b) + half * tol;
            xm = half * (c - b);
            /* bail out if the solution is found to the desired accuracy */
            if ((abs(xm) < tol1) || (fb == zero)) exit = 1;

            if (exit == 0) {  /* proceed with step */

                /* Is bisection necessary ? */
                if ((abs(e) < tol1) || (abs(fa) <= abs(fb))) {
                    d = xm;          /* bisection */
                    e = d;
                } else {
                    if (a == c)  /* use quadratic interp. if are a and c distinct */ {
                        s = fb / fa;          /* linear interpolation */
                        p = two * xm * s;
                        q = one - s;
                    } else {
                        q = fa / fc;          /* quadratic interpolation */
                        r = fb / fc;
                        s = fb / fa;
                        p = s * (two * xm * q * (q - r) - (b - a) * (r - one));
                        q = (q - one) * (r - one) * (s - one);
                    }

                    if (p > zero) q = -q;     /* adjust signs */
                    p = abs(p);

                    /* Is the interpolation acceptable? */

                    if (((two * p) > (three * xm * q - abs(tol1 * q))) ||
                            (p >= abs(half * e * q))) {
                        d = xm;           /* use bisection */
                        e = d;
                    } else {
                        e = d;            /* use previously selected  */
                        d = p / q;        /* interpolation            */
                    }
                }  /* if .. bisection necessary ? */

                /* Complete step */
                a = b;                   /* save old point b as a */
                fa = fb;
                if (abs(d) > tol1)
                    b = b + d;         /* move b to a new point closer to the zero */
                else {                  /* move b by a relatively small amount */
                    if (xm > zero)
                        b = b + abs(tol1);
                    else
                        b = b - abs(tol1);
                }
                fb = F.apply(b);            /* function value at the new point */

                if ((fb * (fc / abs(fc))) <= zero)
                    skip = 1;       /* zero is already between b and c */
                else skip = 0; /* swap a and c to get zero between b and c */

            }  /* if not exit , end of step */
        }  /* while */


        /* all done */
        flag = zflag;

        /* return the abscissa with the minimum absolute value */
        return b;
    }
}
