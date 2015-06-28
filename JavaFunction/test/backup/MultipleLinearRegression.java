
public class MultipleLinearRegression {
	private final int Nx;        // number of X samples
	private final int Ny;		// number of Y samples
	private int Dx;        // number of dependent variables
	private final double[] beta;  // regression coefficients
	private double SSE;         // residual sum of squared
	private double SST;         // total sum of squares
	private double[] y_hat;		// predicted 
	private double y_mean;		// mean of y
	private int DFE;			// degree of freedom of error
	private double s2; 			// square of estimation of residual standard deviation
	private double[] betaSD; // standard deviation of regression coefficients
	private double[] tStats;	// t statistics for regression coefficients 

	public MultipleLinearRegression(double[][] X, double[] y, boolean biasterm) {
		if (X.length != y.length) throw new RuntimeException("dimensions don't agree");
		Ny = y.length;
		Nx = X.length;
		Dx = X[0].length;
		
		if (biasterm) {
			int Dx_new = Dx + 1;
			
			double[][] Xtemp = new double[Nx][Dx_new];
			for ( int i = 0; i < Nx; i++) {
				for (int j = 0; j < Dx; j++) {
					if ( j == 0) {
						Xtemp[i][j] = 1;
					} else {
						Xtemp[i][j] = X[i][j-1];
					}
				}
			}
			X = Xtemp;
			
			DFE = Nx - Dx_new;
			
		} else {
			DFE = Nx - Dx;
		}
		
		
		
		double [][] Xvar = Matrix.multiply(Matrix.transpose(X),X);
		double [][] invXvar = Matrix.invert(Xvar);
        
		beta = Matrix.multiply(Matrix.multiply(invXvar,Matrix.transpose(X)),y);
		
        double sum = 0.0;
        for (int i = 0; i < Ny; i++) {
            sum += y[i];
		}
        y_mean = sum / Ny;

        // total variation to be accounted for
        for (int i = 0; i < Ny; i++) {
            double dev = y[i] - y_mean;
            SST += dev*dev;
        }
		
		y_hat = Matrix.multiply(X,beta);

        // variation not accounted for
        double[] residuals =  Matrix.subtract(y_hat,y);           
        SSE = Matrix.dot(residuals,residuals);
		s2 = SSE / DFE;
		
		// statistical significance test for regression coefficients
		double[] diagVarBeta = Matrix.diag(Matrix.dotmultiply(s2,invXvar));
		betaSD = new double[DFE];
		tStats = new double[DFE];
		// for (int i = 0; i < DFE; i++) {
			// betaSD[i] = Math.sqrt(diagVarBeta[i]);
			// tStats[i] = beta[i]/betaSD[i];
		// }
		
		
		

    }

    public double[] getbeta() {
        return beta;
    }
		
		public double[] getbetaSD() {
			return betaSD;
		}
		
		public double[] gettStats() {
			return tStats;
		}

    public double R2() {
        return 1.0 - SSE/SST;
    }

	// test client
	
    // public static void main(String[] args) {
        // double[][] x = { {  1,  10,  20 },
                         // {  1,  20,  40 },
                         // {  1,  40,  15 },
                         // {  1,  80, 100 },
                         // {  1, 160,  23 },
                         // {  1, 200,  18 } };
        // MultipleLinearRegression regression = new MultipleLinearRegression(x, y, true);
        // double[] y = { 243, 483, 508, 1503, 1764, 2129 };
		

        
    // }
}
