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
    public static void multiplyMatrix(double[][] first_matrix, double[][] second_matrix, double[][] result) {
        int n = first_matrix.length; // Assuming square matrices of size n x n

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < n; k++) {
                    sum += first_matrix[i][k] * second_matrix[k][j];
                }
                result[i][j] = sum;
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
