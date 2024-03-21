package JAVA.SPLINE_SEVAL;

public class Seval {
    private int n;
    private double u;
    private double[] x;
    private double[] y;
    private double[] b;
    private double[] c;
    private double[] d;
    private int last;

    public Seval(int n, double u, double[] x, double[] y, double[] b, double[] c, double[] d, int last) {
        this.n = n;
        this.u = u;
        this.x = x;
        this.y = y;
        this.b = b;
        this.c = c;
        this.d = d;
        this.last = last;
    }

    public double seval() {
        int i, j, k;
        double w;

        i = last;
        if (i >= n-1) {
            i = 0;
        }
        if (i < 0) {
            i = 0;
        }

        if ((x[i] > u) || (x[i+1] < u)) {  /* ---- perform a binary search ---- */
            i = 0;
            j = n;
            do {
                k = (i + j) / 2;         /* split the domain to search */
                if (u < x[k])  j = k;    /* move the upper bound */
                if (u >= x[k]) i = k;    /* move the lower bound */
            }                        /* there are no more segments to search */
            while (j > i+1);
        }
        last = i;

        /* ---- Evaluate the spline ---- */
        w = u - x[i];
        w = y[i] + w * (b[i] + w * (c[i] + w * d[i]));
        return w;
    }
}
