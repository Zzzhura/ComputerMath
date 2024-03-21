package JAVA.DECOMP_SOLVE;

public class Solve {
    private int n;
    private int ndim;
    private double[][] a;
    private double[] b;
    private int[] pivot = new int[8];

    public Solve(int n, int ndim, double[][] a, double[] b, int[] pivot) {
        this.n = n;
        this.ndim = ndim;
        this.a = a;
        this.b = b;
        System.arraycopy(pivot, 0, this.pivot, 0, pivot.length);
    }

    public int solve () {
        int    i, j, k, m;
        double t;

        if (n == 1)
        {
            /* trivial */
            b[0] /= a[0][0];
        }
        else
        {
            /* Forward elimination: apply multipliers. */
            for (k = 0; k < n-1; k ++)
            {
                m = pivot[k];
                t = b[m]; b[m] = b[k]; b[k] = t;
                for (i = k+1; i < n; ++i) {
                    b[i] += a[i][k] * t;
                }
            }

            /* Back substitution. */
            for (k = n-1; k >= 0; --k)
            {
                t = b[k];
                for (j = k+1; j < n; ++j) {
                    t -= a[k][j] * b[j];
                }
                b[k] = t / a[k][k];
            }
        }

        return 0;
    }
}
