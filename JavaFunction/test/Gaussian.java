public class Gaussian {

    // return GaussianPdf(x) = standard Gaussian pdf
    public static double GaussianPdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return GaussianPdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double GaussianPdf(double x, double mu, double sigma) {
        return GaussianPdf((x - mu) / sigma) / sigma;
    }

    // return GaussianCdf(z) = standard Gaussian cdf using Taylor approximation
    public static double GaussianCdf(double z) {
        if (z < -8.0) 
					return 0.0;
        if (z >  8.0) 
					return 1.0;
        double sum = 0.0, term = z;
        for (int i = 3; sum + term != sum; i += 2) {
            sum  = sum + term;
            term = term * z * z / i;
        }
        return 0.5 + sum * GaussianPdf(z);
    }

    // return GaussianCdf(z, mu, sigma) = Gaussian cdf with mean mu and stddev sigma
    public static double GaussianCdf(double z, double mu, double sigma) {
        return GaussianCdf((z - mu) / sigma);
    } 

    // Compute z such that Phi(z) = y via bisection search
    public static double PhiInverse(double y) {
        return PhiInverse(y, .00000001, -8, 8);
    } 

    // bisection search
    private static double PhiInverse(double y, double delta, double lo, double hi) {
        double mid = lo + (hi - lo) / 2;
        if (hi - lo < delta) return mid;
        if (GaussianCdf(mid) > y) 
					return PhiInverse(y, delta, lo, mid);
        else              
					return PhiInverse(y, delta, mid, hi);
    }



    // test client
    // public static void main(String[] args) {
        // double z     = Double.parseDouble(args[0]);
        // double mu    = Double.parseDouble(args[1]);
        // double sigma = Double.parseDouble(args[2]);
        // StdOut.println(Phi(z, mu, sigma));
        // double y = Phi(z);
        // StdOut.println(PhiInverse(y));
    // }

}
