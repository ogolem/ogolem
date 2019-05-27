package contrib.newuoa;


/**
 * Testing routine for the NEWUOA algorithm.
 * @author Jonas Feldt
 * @author Johannes Dieterich
 * @version 2011-11-19
 */
class TestNewUO {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		final double[] x = new double[10];
		final int maxFun = 5000;
		final double rhoEnd = 1E-6;
		int npt;
		for (int n = 2; n < 3; n += 2) {
			npt = 2 * n + 1;
			for (int i = 0; i < n; i++) {
				x[i] = (double) (i + 1) / (double) (n + 1);
			}
			double rhoBeg = 0.2 * x[0];
			System.out.println("Results with N = " + n + " and npt = " + npt);
			final NEWUOAMethod met = new ChebyQuad();
			final NEWUOAOptimizer opt = new NEWUOAOptimizer(rhoBeg, rhoEnd, maxFun);
			opt.doOptimize(n, x, met);			
		}
	}

}


/**
 * An example problem for the NEWUOA optimization.
 * @author Jonas Feldt
 * @author Johannes Dieterich
 * @version 2011-11-19
 */
class ChebyQuad implements NEWUOAMethod {

	@Override
	public double computeObjectiveValue(final int n, double[] point) {
		final double[][] y = new double[10][10];
//		final int n = point.length;
		for (int j = 0; j < n; j++) {
			y[0][j] = 1;
			y[1][j] = 2 * point[j] - 1;			
		}
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < n; j++) {
				y[i + 1][j] = 2.0 * y[1][j] * y[i][j] - y[i - 1][j];
			}
		}
		double f = 0;
		int np = n + 1;
		int iw = 1;
		double sum = 0;
		for (int i = 0; i < np; i++) {
			sum = 0;
			for (int j = 0; j < n; j++) {
				sum += y[i][j];
			}
			sum /= (double) n;
			if (iw > 0) {
				final int t = i + 1;
				sum += 1.0 / (double) (t * t - 2 * t);
			}
			iw = -iw;
			f += sum * sum;
		}
		return f;
	}
}
