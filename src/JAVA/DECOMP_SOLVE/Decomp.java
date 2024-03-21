package JAVA.DECOMP_SOLVE;

import static java.lang.Math.abs;

public class Decomp {
    private double EPSILON = 2.2e-16;

    private int n;
    private int ndim;
    private double[][] a;
    private double cond;
    private int[] pivot;
    private int flag;

    public Decomp(int n, int ndim, double[][] a, double cond, int[] pivot, int flag) {
        this.n = n;
        this.ndim = ndim;
        this.a = a;
        this.cond = cond;
        this.pivot = pivot;
        this.flag = flag;
    }

    public double getCond() {
        return cond;
    }

    public void getMatrix()
    {
        for (var i: a){
            for (int d = 0; d < i.length; ++d)
            {
                System.out.println(i[d]+ ' ');
            }
            System.out.println("\n");
        }
    }
    public int getFlag() {
        return flag;
    }

    public int decomp() {
        double ek, t, pvt, anorm, ynorm, znorm;
        int i, j, k, m;
        double pa, pb;
        double[] work;

        flag = 0;
        work = null;
        if (a == null || pivot == null || n < 1 || ndim < n)
        {
            flag = 2;
            return 0;
        }

        pivot[n-1] = 1;
        if (n == 1)
        {
            /* One element only */
            cond = 1.0;
            if (a[0][0] == 0.0)
            {
                cond = 1.0e+32;  /* singular */
                flag = 3;
                return 0;
            }
            return (0);
        }

        work = new double[n * 8];
        if (work == null)
        {
            flag = 1;
            return 0;
        }

        /* --- compute 1-norm of a --- */

        anorm = 0.0;
        for (j = 0; j < n; ++j)
        {
            t = 0.0;
            for (i = 0; i < n; ++i) {
                t += abs(a[i][j]);
            }
            if (t > anorm) {
                anorm = t;
            }
        }

        /* Apply Gaussian elimination with partial pivoting. */

        for (k = 0; k < n-1; ++k)
        {
   /* Find pivot and label as row m.
      This will be the element with largest magnitude in
      the lower part of the kth column. */
            m = k;
            pvt = abs(a[m][k]);
            for (i = k+1; i < n; ++i)
            {
                t = abs(a[i][k]);
                if ( t > pvt )  {
                    m = i;
                    pvt = t;
                }
            }
            pivot[k] = m;
            pvt = a[m][k];

            if (m != k)
            {
                pivot[n-1] = -pivot[n-1];
                /* Interchange rows m and k for the lower partition. */
                for (j = k; j < n; ++j)
                {
                    t = a[m][j];
                    a[m][j] = a[k][j];
                    a[k][j] = t;
                }
            }
            /* row k is now the pivot row */

            /* Bail out if pivot is too small */
            if (abs(pvt) < anorm * EPSILON)
            {
                /* Singular or nearly singular */
                cond = 1.0e+32;
                flag = 3;
                return 0;
            }

   /* eliminate the lower matrix partition by rows
      and store the multipliers in the k sub-column */
            for (i = k+1; i < n; ++i)
            {
                pa = a[i][k];          /* element to eliminate */
                t = -(pa / pvt);          /* compute multiplier   */
                a[i][k] = t;                     /* store multiplier     */
                for (j = k+1; j < n; ++j)    /* eliminate i th row */
                {
                    if (abs(t) > anorm * EPSILON)
                        a[i][j] += a[k][j] * t;
                }
            }

        }  /* End of Gaussian elimination. */

/* cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
   estimate obtained by one step of inverse iteration for the
   small singular vector. This involves solving two systems
   of equations, (a-transpose)*y = e and a*z = y where e
   is a vector of +1 or -1 chosen to cause growth in y.
   estimate = (1-norm of z)/(1-norm of y)
   Solve (a-transpose)*y = e   */

        for (k = 0; k < n; ++k)
        {
            t = 0.0;
            if (k != 0)
            {
                for (i = 0; i < k; ++i) {
                    t += a[i][k] * work[i];
                }
            }
            if (t < 0.0) {
                ek = -1.0;
            }
            else {
                ek = 1.0;
            }
            pa = a[k][k];
            if (abs(pa) < anorm * EPSILON)
            {
                /* Singular */
                cond = 1.0e+32;
                flag = 3;
                return 0;
            }

            work[k] = -(ek + t) / pa;
        }

        for (k = n-2; k >= 0; --k)
        {
            t = 0.0;
            for (i = k+1; i < n; i++)
                t += a[i][k] * work[i];
      /* we have used work[i] here, however the use of work[k]
	 makes some difference to cond */
            work[k] = t;
            m = pivot[k];
            if (m != k) {
                t = work[m];
                work[m] = work[k];
                work[k] = t;
            }
        }

        ynorm = 0.0;
        for (i = 0; i < n; ++i) ynorm += abs(work[i]);

        /* --- solve a * z = y */
        Solve solve = new Solve(n, ndim, a, work, pivot);
        solve.solve();

        znorm = 0.0;
        for (i = 0; i < n; ++i) {
            znorm += abs(work[i]);
        }

        /* --- estimate condition --- */
        cond = anorm * znorm / ynorm;
        if (cond < 1.0) {
            cond = 1.0;
        }
        if (cond + 1.0 == cond) {
            flag = 3;
        }
        return 0;
    }   /* --- end of function decomp() --- */

}

