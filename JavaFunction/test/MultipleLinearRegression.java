
public class MultipleLinearRegression {
	private final int Nx;        // number of X samples
	private final int Ny;		// number of Y samples
	private int Dx;        // number of dependent variables
	private final double[] beta;  // regression coefficients
	private double SSE;         // residual sum of squared
	private double SST;         // total sum of squares
	private double SSR;					// sum of squares of regression
	private double[] y_hat;		// predicted 
	private double y_mean;		// mean of y
	private int DFE;			// degree of freedom of error
	private int DFR; 			// degree of freedom of regression
	private int DFT; 			// total degree of freedom 
	private double s2; 			// square of estimation of residual standard deviation
	private double[] betaSD; // standard deviation of regression coefficients
	private double[] tStats;	// t statistics for regression coefficients
	private double R2;				// coefficient of determination
	private double R2bar;			// adjusted R2bar
	private double Fstats;		// F statistics to test if all betas are zero
	
	
	

	public MultipleLinearRegression(double[][] X, double[] y, boolean biasterm) {
		if (X.length != y.length) throw new RuntimeException("dimensions don't agree");
		Ny = y.length;
		Nx = X.length;
		Dx = X[0].length;
		
		if (biasterm) {
			Dx = Dx + 1;
			
			double[][] Xtemp = new double[Nx][Dx];
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
			
			DFE = Nx - Dx;
			DFR = Dx - 1;
			DFT = Nx - 1;
			
		} else {
			DFE = Nx - Dx;
			DFR = Dx;
			DFT = Nx;
		}
		
		
		
		double[][] Xvar = Matrix.multiply(Matrix.transpose(X),X);
		double[][] invXvar = Matrix.invert(Xvar);
        
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
				
				for (int i = 0; i < Ny; i++) {
            double dev = y_hat[i] - y_mean;
            SSR += dev*dev;
        }

        // variation not accounted for
        double[] residuals =  Matrix.subtract(y_hat,y);           
        SSE = Matrix.dot(residuals,residuals);
				s2 = SSE / DFE;
		
		// statistical significance test for regression coefficients
		double[] diagVarBeta = Matrix.diag(Matrix.dotmultiply(s2,invXvar));
		betaSD = new double[Dx];
		tStats = new double[Dx];
		for (int i = 0; i < Dx; i++) {
			
			betaSD[i] = Math.sqrt(diagVarBeta[i]);
		
			tStats[i] = beta[i]/betaSD[i];
			
		}
		
		R2 = 1.0 - SSE/SST;
		R2bar = R2 - (1 -R2)*(DFR/DFE);
		
		Fstats = (SSR/DFR)/(SSE/DFE);
		

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

    public double getR2() {
        return R2;
    }
		
		public double getR2bar() {
			return R2bar;
		}
		
		public double gets2() {
			return s2;
		}
		
		public double getFstats() {
			return Fstats;
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