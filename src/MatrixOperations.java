public class MatrixOperations {
    public static void changeRowsColomns(double[][] arr)
    {
        double[][] tmp = new double[8][8];
        for(int i = 0; i <= 7; i++)
        {
            for(int j = 0; j <= 7; j++)
            {
                tmp[j][i] = arr[i][j];
            }
        }
        System.arraycopy(tmp, 0, arr , 0, arr.length);
    }
    public static void multiplyMatrix(double[][] first_matrix, double[][] second_matrix, double[][] result)
    {
        int k = 0;
        int index_y = 0;
        var res = 0;
        for (int c = 0; c < 8; c++)
        {
            for(int i = 0; i < 8; i++)
            {
                for(int j = 0; j < 8; j++)
                {
                    res += first_matrix[c][j] * second_matrix[i][j];
                }
                result[c][i] = res;
                res = 0;
            }
        }
    }

    public static void printMatrix(double[][] R)
    {
        for (int i = 0; i <= 7; i++)
        {
            for(int j = 0; j <= 7; j++)
            {
                System.out.print(R[i][j] + " ");
            }
            System.out.println();
        }
    }
    public static void subtractMatrix(double[][] first_matrix, double[][] second_matrix, double[][] result)
    {
        for(int i = 0; i < 8;i++)
        {
            for(int j = 0; j < 8; j++)
            {
                result[i][j] = first_matrix[i][j] - second_matrix[i][j];
            }
        }
    }

    public static double[] getArrayOfRowsSum(double[][] matrix)
    {
        double [] result = new double[8];
        int i = 0;
            for(int j = 0; j < 8; j++)
            {
                for(int k = 0; k < 8; k++)
                {
                    result[i] += matrix[j][k];
                }
                i++;
            }
        return result;
    }
}
