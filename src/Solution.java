import JAVA.DECOMP_SOLVE.Decomp;
import JAVA.DECOMP_SOLVE.Solve;

import java.util.Arrays;
import java.util.OptionalDouble;

public class Solution {
    private static final double cond = 0.0;
    private  static final int n = 8;
    private  static final int flag = 0;
    private static final int[] ipvt = {0,0,0,0,1,1,0,0};
    public static void doTask()
    {
        double[] p_values = new double[]{1.0, 0.1, 0.01, 0.0001, 0.000001};
        for (double p: p_values)
        {
            double[][] A = {
                    {p + 13, 2, 8, -7, 7, 5, -7, -7},
                    {7, 2, -4, 2, 3, 3, -1, -2},
                    {-7, 2, 1, 3, 6, -6, -3, -4},
                    {-2, -8, -6, -1, 6, 2, 1, -4},
                    {0, 4, -7, 1, 22, 0, -6, -6},
                    {0, -3, -6, 6, 4, 13, 0, 6},
                    {-8, -6, -4, 7, -5, -5, -2, 1},
                    {5, 5, -2, -2, -3, 0, -7, 14}
            };
            double[][] A_copy = {
                    {p + 13, 2, 8, -7, 7, 5, -7, -7},
                    {7, 2, -4, 2, 3, 3, -1, -2},
                    {-7, 2, 1, 3, 6, -6, -3, -4},
                    {-2, -8, -6, -1, 6, 2, 1, -4},
                    {0, 4, -7, 1, 22, 0, -6, -6},
                    {0, -3, -6, 6, 4, 13, 0, 6},
                    {-8, -6, -4, 7, -5, -5, -2, 1},
                    {5, 5, -2, -2, -3, 0, -7, 14}
            };
            double[][] E = {{1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1},
            };
            Decomp decomp = new Decomp(n, n, A, cond, ipvt, flag);
            decomp.decomp();
            for (int i = 0; i <= 7; i++)
            {
                Solve solve = new Solve(n, n, A, E[i], ipvt);
                solve.solve();
            }
            double[][] newArr = {{1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1},
            };
            System.arraycopy(E, 0, newArr, 0, E.length);
            newArr = E;
            MatrixOperations.changeRowsColomns(newArr);
            double[][] R = new double[8][8];
            MatrixOperations.multiplyMatrix(A_copy, newArr, R);
            double[][] Result = new double[8][8];
            double[][] E_1= {{1, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0},
                    {0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 1},
            };
            MatrixOperations.subtractMatrix(R, E_1, Result);
            double[] arrayOfResults = MatrixOperations.getArrayOfRowsSum(Result);
            OptionalDouble result = Arrays.stream(arrayOfResults).max();
            String formattedValue = String.format("%.16f", result.getAsDouble());
            System.out.println("For p = " + p + " RESULT: " + formattedValue);
        }
    }

}
