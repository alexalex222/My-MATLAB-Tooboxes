public class Matrix {

    // return a random m-by-n matrix with values between 0 and 1
    public static double[][] random(int m, int n) {
        double[][] C = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = Math.random();
        return C;
    }

    // return n-by-n identity matrix I
    public static double[][] identity(int n) {
        double[][] I = new double[n][n];
        for (int i = 0; i < n; i++)
            I[i][i] = 1;
        return I;
    }
		
		

    // return x^T y
    public static double dot(double[] x, double[] y) {
        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double sum = 0.0;
        for (int i = 0; i < x.length; i++)
            sum += x[i] * y[i];
        return sum;
    }

    // return C = A^T
    public static double[][] transpose(double[][] A) {
        int m = A.length;
        int n = A[0].length;
        double[][] C = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[j][i] = A[i][j];
        return C;
    }

    // return C = A + B
    public static double[][] add(double[][] A, double[][] B) {
        int m = A.length;
        int n = A[0].length;
        double[][] C = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }
	
	// return C = A - B
    public static double[] subtract(double[] A, double[] B) {
        int m = A.length;
        if (A.length != B.length) throw new RuntimeException("Matrix dimension mismatch");
        double[] C = new double[m];
        for (int i = 0; i < m; i++)
			C[i] = A[i] - B[i];
        return C;
    }

    // return C = A - B
    public static double[][] subtract(double[][] A, double[][] B) {
        int m = A.length;
        int n = A[0].length;
		if ((A.length != B.length) || (A[0].length != B[0].length))  throw new RuntimeException("Matrix dimension mismatch");
        double[][] C = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    // return C = A * B
    public static double[][] multiply(double[][] A, double[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        if (nA != mB) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] C = new double[mA][nB];
        for (int i = 0; i < mA; i++)
            for (int j = 0; j < nB; j++)
                for (int k = 0; k < nA; k++)
                    C[i][j] += (A[i][k] * B[k][j]);
        return C;
    }
		
		// return C = a .* B
    public static double[][] dotmultiply(double a, double[][] B) {
        int mB = B.length;
        int nB = B[0].length;
        double[][] C = new double[mB][nB];
        for (int i = 0; i < mB; i++)
            for (int j = 0; j < nB; j++)
							C[i][j] = a*B[i][j]; 
        return C;
    }
		
		// return C = a .* B
    public static double[][] dotmultiply(double[][] B, double a) {
        int mB = B.length;
        int nB = B[0].length;
        double[][] C = new double[mB][nB];
        for (int i = 0; i < mB; i++)
            for (int j = 0; j < nB; j++)
							C[i][j] = a*B[i][j]; 
        return C;
    }
		
		// return C = A .* B
    public static double[][] dotmultiply(double[][] A, double[][] B) {
				int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
				if ((mA != mB) || (nA != nB)) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] C = new double[mB][nB];
        for (int i = 0; i < mB; i++)
            for (int j = 0; j < nB; j++)
							C[i][j] = A[i][j]*B[i][j]; 
        return C;
    }
	
	public static double[] diag (double[][] A) {
		int rows = A.length;
		int columns = A[0].length;
		int n = (int) Math.min(rows, columns);
		double [] C = new double [n];
		for (int i = 0; i < n; i++) {
			C[i] = A[i][i];
		}
		return C;
	}

    // matrix-vector multiplication (y = A * x)
    public static double[] multiply(double[][] A, double[] x) {
        int m = A.length;
        int n = A[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i] += (A[i][j] * x[j]);
        return y;
    }


    // vector-matrix multiplication (y = x^T A)
    public static double[] multiply(double[] x, double[][] A) {
        int m = A.length;
        int n = A[0].length;
        if (x.length != m) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j] += (A[i][j] * x[i]);
        return y;
    }
	
	public static double[][] invert(double[][] A){
		int m = A.length;
        int n = A[0].length;
		if (m != n) throw new RuntimeException("Not a square matrix");
		double[][] X = new double[n][n];
		double[][] B = new double[n][n];
		int[] index = new int[n];
		for (int i=0; i<n; ++i) {
			B[i][i] = 1;
		}
		
		// Transform the matrix into a upper triangle Gaussian partial pivot
		
		gaussianPP(A, index);
		
		// Update the matrix B[i][j] with the ratios stored 
		for (int i = 0; i < n-1; ++i) {
			for(int j= i+1; j < n; ++j) {
				for(int k = 0; k < n; ++k) {
					B[index[j]][k] -= A [index[j]][i]*B[index[i]][k];
				}
			}
		}
		
		for (int i = 0; i<n; ++i){
			X[n-1][i] = B[index[n-1]][i]/A[index[n-1]][n-1];
			for (int j=n-2; j>=0; --j) {
				X[j][i] = B[index[j]][i];
				for (int k=j+1; k<n; ++k) {
					X[j][i] -= A[index[j]][k]*X[k][i];
				}
				X[j][i] /= A[index[j]][j];
			}
		}
			
		return X;					
		
	}
	
	// Method to carry out the Gaussian partial pivot elimination. Here index[] stores pivoting order
	
	public static void gaussianPP(double[][] A, int[] index){
		int n = index.length;
		double[] c = new double[n];
		
		// initialize the index
		for (int i = 0; i < n; ++i){
			index[i] = i;
		}
		
		// find the scaling factor, one from each row
		for (int i=0; i<n; ++i) {
			double c1 = 0;
			for (int j = 0 ; j<n; ++j) {
				double c0 = Math.abs(A[i][j]);
				if (c0 > c1) c1 = c0;
			}
			c[i] = c1;				
		}
		
		// search the pivoting element from each column
		int k =0;
		for (int j = 0; j<n-1; ++j) {
			double pi1 = 0;
			for (int i = j; i<n; ++i) {
				double pi0 = Math.abs(A[index[i]][j]);
				pi0/=c[index[i]];
				if (pi0 > pi1) {
				pi1 = pi0;
				k = i;
				}
			}
			// interchange rows according to the pivoting order
			int itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;
			for (int i=j+1; i<n; ++i) {
				double pj = A[index[i]][j]/A[index[j]][j];
				A[index[i]][j] = pj;
				
				for (int l = j+1; l<n; ++l) {
					A[index[i]][l] -= pj*A[index[j]][l];
				}
			}
		}		
	}
	
    
}