package contrib.newuoa;

import java.io.Serializable;
import org.ogolem.core.FixedValues;
import org.ogolem.helpers.Tuple;

/**
 * The actual NEWUAO optimization engine.
 * @author Jonas Feldt
 * @author Johannes Dieterich
 * @version 2012-06-24
 */
public class NEWUOAOptimizer implements Serializable {
	
    private static final long serialVersionUID = (long) 20111119;
    private static final double TWO_PI = 2.0*Math.PI;
	
    final double rhoBeg;
    final double rhoEnd;
    final int maxFun;
	
	public NEWUOAOptimizer(final double rhoBeg, final double rhoEnd,
			final int maxFun) {
		this.rhoBeg = rhoBeg;
		this.rhoEnd = rhoEnd;
		this.maxFun = maxFun;
	}
        
        public NEWUOAOptimizer(final NEWUOAOptimizer orig){
            this.rhoBeg = orig.rhoBeg;
            this.rhoEnd = orig.rhoEnd;
            this.maxFun = orig.maxFun;
        }

	/**
	 * Do an optimization using NEWUAO.
	 * @param guess
	 * @param met
	 * @return a tupel containing the function value and the optimized solution vector
	 */
	public Tuple<Double, double[]> doOptimize(final int n, final double[] guess,
			final NEWUOAMethod met) {
//		final int n = guess.length;
		final int npt = 2 * n + 1;
		// DEBUG
		final short iPrint = 0;
                double fit = FixedValues.NONCONVERGEDENERGY;
		try {
			fit = newuoa(n, npt, guess, rhoBeg, rhoEnd, iPrint, maxFun, met);
		} catch (Exception e) {
			e.printStackTrace(System.err);
                        return new Tuple<>(FixedValues.NONCONVERGEDENERGY, guess);
		}
		return new Tuple<>(fit,guess);
	}
	
	/**
     This subroutine seeks the least value of a function of many variables,
     by a trust region method that forms quadratic models by interpolation.
     There can be some freedom in the interpolation conditions, which is
     taken up by minimizing the Frobenius norm of the change to the second
     derivative of the quadratic model, beginning with a zero matrix. The
     arguments of the subroutine are as follows.

     @param N must be set to the number of variables and must be at least two.
     @param NPT is the number of interpolation conditions. Its value must be in the
       interval [N+2,(N+1)(N+2)/2].
     @param x Initial values of the variables must be set in X(1),X(2),...,X(N). They
       will be changed to the values that give the least calculated F.
     @param rhoBeg
     @param rhoEnd RHOBEG and RHOEND must be set to the initial and final values of a trust
       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
       RHOBEG should be about one tenth of the greatest expected change to a
       variable, and RHOEND should indicate the accuracy that is required in
       the final values of the variables.
     @param iPrint The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
       amount of printing. Specifically, there is no output if IPRINT=0 and
       there is output only at the return if IPRINT=1. Otherwise, each new
       value of RHO is printed, with the best vector of variables so far and
       the corresponding value of the objective function. Further, each new
       value of F with its variables are output if IPRINT=3.
     @param MAXFUN must be set to an upper bound on the number of calls of CALFUN.
     
     
     The array W will be used for working space. Its length must be at least
     (NPT+13)*(NPT+N)+3*N*(N+3)/2.

     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
     the value of the objective function for the variables X(1),X(2),...,X(N).

     @param w Partition the working space array, so that different parts of it can be
     treated separately by the subroutine that performs the main calculation.
     *
	 * @throws Exception Return from NEWUOA because NPT is not in the required interval.
	 */
	private double newuoa(final int n, final int npt, final double[] x,
			final double rhoBeg, final double rhoEnd, final short iPrint,
			final int maxFun, final NEWUOAMethod met) throws Exception {
		
		final int np = 	n + 1;
		final int nptm = npt - np;
		
		if (npt < (n + 2) || npt > ((n + 2) * np) / 2) {
			throw new Exception("Return from NEWUOA because NPT is not in the required interval:" + npt);
		}
		
		final int nDim = npt + n;

	//	The above settings provide a partition of W for subroutine NEWUOB.
	//	The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
	//	W plus the space that is needed by the last array of NEWUOB.
		
		final double[] ixb = new double[n];
		final double[] ixo = new double[n];
		final double[] ixn = new double[n];
		final double[][] ixp = new double[npt][n];
		final double[] ifv = new double[npt];
		final double[] igq = new double[n];
		final double[] ihq = new double[(n * np) / 2];
		final double[] ipq = new double[npt];
		final double[][] bMat = new double[nDim][n];
		final double[][] zMat = new double[npt][nptm];
		final double[] id = new double[n];
		final double[] ivl = new double[nDim];
		final ScopedPtr w = new ScopedPtr(new double[2 * npt], 0);
		final double[][] wVec = new double[nDim][5];
		final double[][] prod = new double[nDim][5];

		return newuob(n, npt, x, rhoBeg, rhoEnd, iPrint, maxFun, ixb, ixo, ixn,
				ixp, ifv, igq, ihq, ipq, bMat, zMat, nDim, id, ivl, w, wVec, prod, met);
		
	}
	
	
	/**
	 * The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
       to the corresponding arguments in SUBROUTINE NEWUOA.
     @param XBASE will hold a shift of origin that should reduce the contributions
       from rounding errors to values of the model and Lagrange functions.
     @param XOPT will be set to the displacement from XBASE of the vector of
       variables that provides the least calculated F so far.
     @param XNEW will be set to the displacement from XBASE of the vector of
       variables for the current calculation of F.
     @param XPT will contain the interpolation point coordinates relative to XBASE.
     @param FVAL will hold the values of F at the interpolation points.
     @param GQ will hold the gradient of the quadratic model at XBASE.
     @param HQ will hold the explicit second derivatives of the quadratic model.
     @param PQ will contain the parameters of the implicit second derivatives of
       the quadratic model.
     @param BMAT will hold the last N columns of H.
     @param ZMAT will hold the factorization of the leading NPT by NPT submatrix of
       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
       the elements of DZ are plus or minus one, as specified by IDZ.
     @param NDIM is the first dimension of BMAT and has the value NPT+N.
     @param d is reserved for trial steps from XOPT.
     @param VLAG will contain the values of the Lagrange functions at a new point X.
       They are part of a product that requires VLAG to be of length NDIM.
     @param w The array W will be used for working space. Its length must be at least
       10*NDIM = 10*(NPT+N).
	 */
        @SuppressWarnings("fallthrough")
	private double newuob(final int n, final int npt, final double[] x,
			final double rhoBeg, final double rhoEnd, final short iPrint,
			final int maxFun, final double[] xBase, final double[] xOpt,
			final double[] xNew, final double[][] xpt, final double[] fVal,
			final double[] gq, final double[] hq, final double[] pq,
			final double[][] bMat, final double[][] zMat, final int nDim,
			final double[] d, final double[] vLag, final ScopedPtr w,
			final double[][] wVec, final double[][] prod, final NEWUOAMethod met) {
		

		// Set some constants.
		
		final int np = n;
		final int nh = (n * np) / 2 + 1;
		final int nptm = npt - np - 1;
		final int nfTest = Math.max(maxFun, 1);
                System.arraycopy(x, 0, xBase, 0, n);
		
		double rhoSq = rhoBeg * rhoBeg;
		final double recip = 1.0 / rhoSq;
		final double reciq = Math.sqrt(0.5) / rhoSq;
		int nf = -1; // -1!
		
		double f = 0;
		double fBeg = 0;
		double fOpt = 0;
		int kOpt = 0;
		int nfm = 0;
		int nfmm = 0;
		int nfSav = 0;
		int ipt = 0;
		int jpt = 0;
		double xipt = 0;
		double xjpt = 0;
		double rho = 0;
		double[] crvMin = {0};
		double delta = 0;
		double ratio = 0;
		double diff = 0;
		double diffA = 0;
		double diffB = 0;
		double diffC = 0;
		double xOptSq = 0;
		double dSq = 0;
		double tempQ = 0;
		double sum = 0;
		int idz = 0;
		double temp = 0;
		int ip = 0;
		int ih = 0;
		double dStep = 0;
		double[] alpha = new double[1];
		int kNew = -1;
		double dNorm = 0;
		double beta = 0;
		int kSave = -1;
		int iTest = 0;
		double fSave = 0;
		double vQuad = 0;
		
		int state = 50;		
		for (;;) {
			switch (state) {
		
			case 50:
				nfm = nf;
				nfmm = nf - n;
				nf++;
				if (nfm < 2 * n) {
					if (nfm > -1 && nfm < n) {
						xpt[nf][nfm] = rhoBeg;
					} else if (nfm >= n) {
						xpt[nf][nfmm] = -rhoBeg;
					}
				} else {
					int iTemp = (nfmm - 1) / (n - 1);
					jpt = nfm - iTemp * (n - 1) - (n - 1);
					ipt = jpt + iTemp;
					if (ipt > n) {
						iTemp = jpt;
						jpt = ipt - (n - 1);
						ipt = iTemp;
					}
					xipt = rhoBeg;
					if (fVal[ipt + np] < fVal[ipt + 1]) {
						xipt = -xipt;
					}
					xjpt = rhoBeg;
					if (fVal[jpt + np] < fVal[jpt + 1]) {
						xjpt = -xjpt;
						xpt[nf][ipt] = xipt;
						xpt[nf][jpt] = xjpt;
					}
				}
			
//				Calculate the next value of F, label 70 being reached immediately
//				after this calculation. The least function value so far and its index
//				are required.
				
				for (int j = 0; j < n; j++) {
					x[j] = xpt[nf][j] + xBase[j];
				}
				state = 310;
				continue;
				
			case 70:
				fVal[nf] = f;
				if (nf == 0) { // 0!
					fBeg = f;
					fOpt = f;
					kOpt = 0;
				} else if (f < fOpt) {
					fOpt = f;
					kOpt = nf;
				}
				
//				Set the nonzero initial elements of BMAT and the quadratic model in
//				the cases when NF is at most 2*N+1.
				
				if (nfm < 2 * n) {
					if (nfm > -1 && nfm < n) {
						gq[nfm] = (f - fBeg) / rhoBeg;
						if (npt < nf + n - 1) {
							bMat[0][nfm] = -1.0 / rhoBeg;
							bMat[nf][nfm] = 1.0 / rhoBeg;
							bMat[npt + nfm][nfm] = -0.5 / rhoSq;
						}
					} else if (nfm >= n) {						
						bMat[nf - n][nfmm] = 0.5 / rhoBeg;
						bMat[nf][nfmm] = -0.5 / rhoBeg;
						zMat[0][nfmm] = -2.0 * reciq;
						zMat[nf - n][nfmm] = reciq;
						zMat[nf][nfmm] = reciq;
						ih = (((nfmm + 1) * (nfmm + 2)) / 2) - 1;
						temp = (fBeg - f) / rhoBeg;
						hq[ih] = (gq[nfmm] - temp) / rhoBeg;
						gq[nfmm] = 0.5 * (gq[nfmm] + temp);
					}
					
//					Set the off-diagonal second derivatives of the Lagrange functions and
//					the initial quadratic model.

				} else {
					ih = (ipt * (ipt - 1)) / 2 + jpt - 1;
					if (xipt < 0) {
						ipt += (n - 1);
					}
					if (xjpt < 0) {
						xjpt += (n - 1);
					}
					zMat[0][nfmm] = recip;
					zMat[nf][nfmm] = recip;
					zMat[ipt + 1][nfmm] = -recip;
					zMat[jpt + 1][nfmm] = -recip;
					hq[ih] = (fBeg - fVal[ipt + 1] - fVal[jpt + 1] + f) / (xipt * xjpt);
				}
				if (nf < npt - 1) {
					state = 50;
					continue;
				}				
				
				
//				Begin the iterative procedure, because the initial model is complete.
				
				rho = rhoBeg;
				delta = rho;
				idz = 1;
				diffA = 0;
				diffB = 0;
				iTest = 0;
				xOptSq = 0;
				for (int i = 0; i < n; i++) {
					xOpt[i] = xpt[kOpt][i];
					xOptSq += xOpt[i] * xOpt[i];
				}
				
			case 90:
				nfSav = nf;

//				Generate the next trust region step and test its length. Set KNEW
//				to -1 if the purpose of the next F will be to improve the model.

			case 100:
				kNew = -1;
				trsapp(n, npt, xOpt, xpt, gq, hq, pq, delta, d, w, w.ptr(np),
						w.ptr(np + n), w.ptr(np + 2 * n), crvMin);				
				dSq = 0;
				for (int i = 0; i < n; i++) {
					dSq += d[i] * d[i];
				}
				dNorm = Math.min(delta, Math.sqrt(dSq)); 
				if (dNorm < 0.5 * rho) {
					kNew = -2;
					delta *= 0.1;
					ratio = -1;
					if (delta <= 1.5 * rho) {
						delta = rho;
					}
					if (nf <= nfSav + 2) {
						state = 460;
						continue;
					}
					temp = 0.125 * crvMin[0] * rho * rho;
					if (temp <= (Math.max(diffC, Math.max(diffA, diffB)))) {
						state = 460;
						continue;
					}
					state = 490;
					continue;
				}
				
//				Shift XBASE if XOPT may be too far from XBASE. First make the changes
//				to BMAT that do not depend on ZMAT.
				
			case 120:
				if (dSq <= 1E-3 * xOptSq) {
					tempQ = 0.25 * xOptSq;
					for (int k = 0; k < npt; k++) {
						sum = 0;
						for (int i = 0; i < n; i++) {
							sum += xpt[k][i] * xOpt[i];
						}
						temp = pq[k] * sum;
						sum -= 0.5 * xOptSq;
						w.set(npt + k, sum);
						for (int i = 0; i < n; i++) {
							gq[i] += temp * xpt[k][i];
							xpt[k][i] -= 0.5 * xOpt[i];
							vLag[i] = bMat[k][i];
							w.set(i, sum * xpt[k][i] + tempQ * xOpt[i]);
							ip = npt + i;
							for (int j = 0; j < i + 1; j++) {
								bMat[ip][j] += vLag[i] * w.get(j) + w.get(i) * vLag[j];
							}
						}
					}

//				Then the revisions of BMAT that depend on ZMAT are calculated.				

					for (int k = 0; k < nptm; k++) {
						double sumZ = 0;
						for (int i = 0; i < npt; i++) {
							sumZ += zMat[i][k];
							w.set(i, w.get(i + npt) * zMat[i][k]);
						}
						for (int j = 0; j < n; j++) {
							sum = tempQ * sumZ * xOpt[j];
							for (int i = 0; i < npt; i++) {
								sum += w.get(i) * xpt[i][j];
							}
							vLag[j] = sum;
							if (k < idz - 1) {
								sum = -sum;
							}
							for (int i = 0; i < npt; i++) {
								bMat[i][j] += sum * zMat[i][k];
							}
						}
						for (int i = 0; i < n; i++) {
							ip = i + npt;
							temp = vLag[i];
							if (k < idz - 1) {
								temp = -temp;
							}
							for (int j = 0; j <= i; j++) {
								bMat[ip][j] += temp * vLag[j];
							}
						}
					}

//				The following instructions complete the shift of XBASE, including
//				the changes to the parameters of the quadratic model.
					
					ih = -1;
					for (int j = 0; j < n; j++) {
						w.set(j, 0);
						for (int k = 0; k < npt; k++) {
							w.set(j, w.get(j) + pq[k] * xpt[k][j]);
							xpt[k][j] -= 0.5 * xOpt[j];
						}
						for (int i = 0; i <= j; i++) {							
							ih++;
							if (i < j) {
								gq[j] += hq[ih] * xOpt[i];
							}
							gq[i] += hq[ih] * xOpt[j];
							hq[ih] += w.get(i) * xOpt[j] + xOpt[i] * w.get(j);
							bMat[npt + i][j] = bMat[npt + j][i];
						}
					}
					for (int j = 0; j < n; j++) {
						xBase[j] += xOpt[j];
						xOpt[j] = 0;						
					}
					xOptSq = 0;
				}
				
//				Pick the model step if KNEW is positive. A different choice of D
//				may be made later, if the choice of D by BIGLAG causes substantial
//				cancellation in DENOM.

				if (kNew > -1) {
					ScopedPtr pointer = new ScopedPtr(vLag, 0);
					biglag(n, npt, xOpt, xpt, bMat, zMat, idz, nDim, kNew, dStep, d,
							alpha, pointer, pointer.ptr(npt), w, w.ptr(np), w.ptr(np + n));
				}

//				Calculate VLAG and BETA for the current choice of D. The first NPT
//				components of W_check will be held in W.
				
				double sumA;
				double sumB;
				for (int k = 0; k < npt; k++) {
					sumA = 0;
					sumB = 0;
					sum = 0;
					for (int j = 0; j < n; j++) {
						sumA += xpt[k][j] * d[j];
						sumB += xpt[k][j] * xOpt[j];
						sum += bMat[k][j] * d[j];
					}
					w.set(k, sumA * (0.5 * sumA + sumB));
					vLag[k] = sum;
				}				
				beta = 0;
				for (int k = 0; k < nptm; k++) {
					sum = 0;
					for (int i = 0; i < npt; i++) {
						sum += zMat[i][k] * w.get(i);
					}
					if (k < idz - 1) {
						beta += sum * sum;
						sum = -sum;
					} else {
						beta -= sum * sum;						
					}
					for (int i = 0; i < npt; i++) {
						vLag[i] += sum * zMat[i][k];
					}
				}
				double bSum = 0;
				double dx = 0;
				for (int j = 0; j < n; j++) {
					sum = 0;
					for (int i = 0; i < npt; i++) {
						sum += w.get(i) * bMat[i][j];
					}
					bSum += sum * d[j];
					int jp = npt + j;
					for (int k = 0; k < n; k++) {
						sum += bMat[jp][k] * d[k];
					}
					vLag[jp] = sum;
					bSum += sum * d[j];
					dx += d[j] * xOpt[j];
				}
				beta = dx * dx + dSq * (xOptSq + 2 * dx + 0.5 * dSq) + beta - bSum;
				vLag[kOpt] += 1;
				
//				If KNEW is positive and if the cancellation in DENOM is unacceptable,
//				then BIGDEN calculates an alternative model step, XNEW being used for
//				working space.
				
				if (kNew > -1) {
					temp = 1.0 + alpha[0] * beta / vLag[kNew] * vLag[kNew];
					if (Math.abs(temp) <= 0.8) {
						bigDen(n, npt, xOpt, xpt, bMat, zMat, idz, nDim, kOpt,
								kNew, d, w, vLag, beta, xNew, wVec, prod);
					}
				}

//				Calculate the next value of the objective function.
				
			case 290:
				for (int i = 0; i < n; i++) {
					xNew[i] = xOpt[i] + d[i];
					x[i] = xBase[i] + xNew[i];
				}
				nf++;
								
			case 310:
				if (nf > nfTest) {
					nf--;
					if (iPrint > 0) {
//						System.out.println("Return from NEWUOA because CALFUN has been called MAXFUN times.");
						state = 530;
						continue;
					}
				}
				f = met.computeObjectiveValue(n, x);
				if (iPrint == 3) {
//					System.out.println("Function number " + nf + "   f = "
//							+ f + "  The corresponding X is: ");
//					for (int i = 0; i < n; i++) {
//						System.out.print(x[i] + " ");
//					}
//					System.out.println();
				}
				if (nf < npt) {
					state = 70;
					continue;
				}
				if (kNew == -2) {
					state = 530;
					continue;
				}
				
//				Use the quadratic model to predict the change in F due to the step D,
//				and set DIFF to the error of this prediction.
				
				vQuad = 0;
				ih = -1;
				for (int j = 0; j < n; j++) {
					vQuad += d[j] * gq[j];
					for (int i = 0; i <= j; i++) {
						ih++;
						temp = d[i] * xNew[j] + d[j] * xOpt[i];
						if (i == j) {
							temp *= 0.5;
						}
						vQuad += temp * hq[ih];
					}
				}
				for (int k = 0; k < npt; k++) {
					vQuad += pq[k] * w.get(k);
				}
				diff = f - fOpt - vQuad;
				diffC = diffB;
				diffB = diffA;
				diffA = Math.abs(diff);
				if (dNorm > rho) {
					nfSav = nf;
				}
				
//				Update FOPT and XOPT if the new F is the least value of the objective
//				function so far. The branch when KNEW is positive occurs if D is not
//				a trust region step.
				
				fSave = fOpt;
				if (f < fOpt) {
					fOpt = f;
					xOptSq = 0;
					for (int i = 0; i < n; i++) {
						xOpt[i] = xNew[i];
						xOptSq += xOpt[i] * xOpt[i];
					}
				}
				kSave = kNew;
				if (kNew > -1) {
					state = 410;
					continue;
				}

//				Pick the next value of DELTA after a trust region step.
				
				if (vQuad >= 0) {
					if (iPrint > 0) {
//						System.out.println("Return from NEWUOA because a trust region step has failed to reduce Q.");
					}
					state = 530;
					continue;
				}
				ratio = (f - fSave) / vQuad;
				if (ratio <= 0.1) {
					delta = 0.5 * dNorm;
				} else if (ratio <= 0.7) {
					delta = Math.max(0.5 * delta, dNorm);
				} else {
					delta = Math.max(0.5 * delta, 2 * dNorm);
				}
				if (delta <= 1.5 * rho) {
					delta = rho;
				}

//				Set KNEW to the index of the next interpolation point to be deleted.
				
				final double t = Math.max(0.1 * delta, rho);
				rhoSq = t * t;
				int kTemp = -1;
				double detrat = 0;
				if (f >= fSave) {
					kTemp = kOpt;
					detrat = 1;
				}
				for (int k = 0; k < npt; k++) {
					double hDiag = 0;
					for (int j = 0; j < nptm; j++) {
						temp = 1;
						if (j < idz - 1) {
							temp = -1;
						}
						hDiag += temp * zMat[k][j] * zMat[k][j];
					}
					temp = Math.abs(beta * hDiag + vLag[k] * vLag[k]);
					double distSq = 0;
					for (int j = 0; j < n; j++) {
						final double t1 = xpt[k][j] - xOpt[j];
						distSq += t1 * t1;
					}					
					if (distSq > rhoSq) {
						final double t1 = distSq / rhoSq;
						temp *= t1 * t1 * t1;
					}
					if (temp > detrat && k != kTemp) {
						detrat = temp;
						kNew = k;
					}					
				}
				if (kNew == -1) {
					state = 460;
					continue;
				}
				
//				Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
//				can be moved. Begin the updating of the quadratic model, starting
//				with the explicit second derivative term.

			case 410:
				update(n, npt, bMat, zMat, idz, nDim, vLag, beta, kNew, w);
				fVal[kNew] = f;
				ih = -1;
				for (int i = 0; i < n; i++) {
					temp = pq[kNew] * xpt[kNew][i];
					for (int j = 0; j < i + 1; j++) {
						ih++;
						hq[ih] += temp * xpt[kNew][j];
					}
				}
				pq[kNew] = 0;
				
//				Update the other second derivative parameters, and then the gradient
//				vector of the model. Also include the new interpolation point.

				for (int j = 0; j < nptm; j++) {
					temp = diff * zMat[kNew][j];
					if (j < idz - 1) {
						temp = -temp;
					}
					for (int k = 0; k < npt; k++) {
						pq[k] += temp * zMat[k][j];
					}
				}
				double gqSq = 0;
				for (int i = 0; i < n; i++) {
					gq[i] += diff * bMat[kNew][i];
					gqSq += gq[i] * gq[i];
					xpt[kNew][i] = xNew[i];
				}
				
//				If a trust region step makes a small change to the objective function,
//				then calculate the gradient of the least Frobenius norm interpolant at
//				XBASE, and store it in W, using VLAG for a vector of right hand sides.

				if (kSave == -1 && delta == rho) {
					if (Math.abs(ratio) > 1E-2) {
						iTest = 0;
					} else {
						for (int k = 0; k < npt; k++) {
							vLag[k] = fVal[k] - fVal[kOpt];
						}
						double giSq = 0;
						for (int i = 0; i < n; i++) {
							sum = 0;
							for (int k = 0; k < npt; k++) {
								sum += bMat[k][i] * vLag[k];
							}
							giSq += sum * sum;
							w.set(i, sum);
						}
						
//						Test whether to replace the new quadratic model by the least Frobenius
//						norm interpolant, making the replacement if the test is satisfied.

						iTest++;
						if (gqSq > 1E2 * giSq) {
							iTest = 0;
						}
						if (iTest > 2) {
							for (int i = 0; i < n; i++) {
								gq[i] = w.get(i);
							}
							for (ih = 0; ih < nh; ih++) {
								hq[ih] = 0;
							}
							for (int j = 0; j < nptm; j++) {
								w.set(j, 0);
								for (int k = 0; k < npt; k++) {
									w.set(j, w.get(j) * vLag[k] * zMat[k][j]);
								}
								if (j < idz - 1) {
									w.set(j, -w.get(j));
								}
							}
							for (int k = 0; k < npt; k++) {
								pq[k] = 0;
								for (int j = 0; j < nptm; j++) {
									pq[k] += zMat[k][j] * w.get(j);
								}
							}
							iTest = 0;
						}						
					}
				}
				if (f < fSave) {
					kOpt = kNew;
				}
				
//				If a trust region step has provided a sufficient decrease in F, then
//				branch for another trust region calculation. The case KSAVE>0 occurs
//				when the new function value was calculated by a model step.

				if (f <= fSave + 0.1 * vQuad) {
					state = 100;
					continue;
				}
				if (kSave > -1) {
					state = 100;
					continue;
				}
				
//				Alternatively, find out if the interpolation points are close enough
//				to the best point so far.

				kNew = -1;
				
			case 460:
				double distSq = 4.0 * delta * delta;
				for (int k = 0; k < npt; k++) {
					sum = 0;
					for (int j = 0; j < n; j++) {
						final double t1 = xpt[k][j] - xOpt[j];
						sum += t1 * t1;
					}
					if (sum > distSq) {
						kNew = k;
						distSq = sum;
					}
				}
				
//				If KNEW is positive, then set DSTEP, and branch back for the next
//				iteration, which will generate a "model step".

				if (kNew > -1) {
					dStep = Math.max(Math.min(0.1 * Math.sqrt(distSq), 0.5 * delta), rho);
					dSq = dStep * dStep;
					state = 120;
					continue;
				}
				if (ratio > 0) {
					state = 100;
					continue;
				}
				if (Math.max(delta, dNorm) > rho) {
					state = 100;
					continue;
				}
				
//				The calculations with the current value of RHO are complete. Pick the
//				next values of RHO and DELTA.

			case 490:
				if (rho > rhoEnd) {
					delta = 0.5 * rho;
					ratio = rho / rhoEnd;
					if (ratio <= 16) {
						rho = rhoEnd;
					} else if (ratio <= 250) {
						rho = Math.sqrt(ratio) * rhoEnd;
					} else {
						rho *= 0.1;
					}
					delta = Math.max(delta, rho);
					if (iPrint > 1) {
//						System.out.println("New RHO = " + rho + " Number of function values = " + nf);
//						System.out.print("Least value of F = " + fOpt + " The corresponding X is: ");
//						for (int i = 0; i < n; i++) {
//							System.out.print((xBase[i] + xOpt[i]) + " ");
//						}
//						System.out.println();						
					}
					state = 90;
					continue;
				}
				
				
//				Return from the calculation, after another Newton-Raphson step, if
//				it is too short to have been tried before.
				
				if (kNew == -2) {
					state = 290;
					continue;
				}
				
			case 530:
				if (fOpt <= f) {
					for (int i = 0; i < n; i++) {
						x[i] = xBase[i] + xOpt[i];
					}
					f = fOpt;
				}
				if (iPrint > 0) {
//					System.out.println("At the return from NEWUOA Number of function values = " + nf);
//					System.out.print("Least value of F = " + f + " The corresponding X is: ");
//					for (int i = 0; i < n; i++) {
//						System.out.print(x[i] + " ");
//					}
//					System.out.println();
				}
				return f;
			}
		}
	}
	
	
	/**
	 * @param N is the number of variables of a quadratic objective function, Q say.
	   @param The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
	       in order to define the current quadratic model Q.
	   @param DELTA is the trust region radius, and has to be positive.
	   @param STEP will be set to the calculated trial step.
	   @param The arrays D, G, HD and HS will be used for working space.
	   @param CRVMIN will be set to the least curvature of H along the conjugate
	       directions that occur, except that it is set to zero if STEP goes
	       all the way to the trust region boundary.
	
	     The calculation of STEP begins with the truncated conjugate gradient
	     method. If the boundary of the trust region is reached, then further
	     changes to STEP may be made, each one being in the 2D space spanned
	     by the current STEP and the corresponding gradient of Q. Thus STEP
	     should provide a substantial reduction to Q within the trust region.
	     
	     JF: d in this function is not equal to d in newuob! d from
	     newuob is called step here.   
	 */
        @SuppressWarnings("fallthrough")
	private void trsapp(final int n, final int npt, final double[] xOpt,
			final double[][] xPt, final double[] gq, final double[] hq,
			final double[] pq, final double delta, final double[] step,
			final ScopedPtr d, final ScopedPtr g, final ScopedPtr hd,
			final ScopedPtr hs, double[] crvMin) {

//	
//	     Initialization, which includes setting HD to H times XOPT.

		final double delSq = delta * delta; 
		int iterC = 0;
		final int iterMax = n;
		int iterSw = iterMax;
		for (int i = 0; i < n; i++) {
			d.set(i, xOpt[i]);
		}		
		int state = 170;
		
		double qRed = 0;
		double qAdd;
		double dd = 0;
		double ds = 0;
		double ss = 0;
		double gg = 0;
		double ggBeg = 0;
		double bStep = 0;
		double dhd;
		double temp;
		double shs = 0;
		double sg = 0;
		double tempA = 0;
		double tempB = 0;
		
		for (;;) {
			switch (state) {
			
//			Prepare for the first line search.
			
			case 20:
				qRed = 0;
				dd = 0;
				for (int i = 0; i < n; i++) {
					step[i] = 0;
					hs.set(i, 0);
					g.set(i, gq[i] + hd.get(i));
					d.set(i, -g.get(i));
					final double t = d.get(i);
					dd += t * t;
				}
				crvMin[0] = 0;
				if (dd == 0) {
					state = 160;
					continue;
				}
				ds = 0;
				ss = 0;
				gg = dd;
				ggBeg = gg;

//				Calculate the step to the trust region boundary and the product HD.

			case 40:
				iterC++;
				temp = delSq - ss; 
				bStep = temp / (ds + Math.sqrt(ds * ds + dd * temp));
				state = 170;
				continue;
				
			case 50:
				dhd = 0;
				for (int j = 0; j < n; j++) {
					dhd += d.get(j) * hd.get(j);
				}

//				Update CRVMIN and set the step-length ALPHA.
				
				double alpha = bStep;
				if (dhd > 0) {
					temp = dhd / dd;
					if (iterC == 1) {
						crvMin[0] = temp;
					}
					crvMin[0] = Math.min(crvMin[0], temp);
					alpha = Math.min(alpha, gg / dhd);
				}
				qAdd = alpha * (gg - 0.5 * alpha * dhd);
				qRed += qAdd;

//				Update STEP and HS.				
				
				final double ggSav = gg;
				gg = 0;
				for (int i = 0; i < n; i++) {
					step[i] += alpha * d.get(i);
					hs.set(i, hs.get(i) + alpha * hd.get(i));
					final double t = g.get(i) + hs.get(i);
					gg += t * t;
				}

//				Begin another conjugate direction iteration if required.
				
				if (alpha < bStep) {
					if (qAdd <= 0.01 * qRed) {
						state = 160;
						continue;
					}
					if (gg <= 1.0E-4 * ggBeg) {
						state = 160;
						continue;
					}
					if (iterC == iterMax) {
						state = 160;
						continue;
					}
					
					temp = gg / ggSav;
					dd = 0;
					ds = 0;
					ss = 0;
					for (int i = 0; i < n; i++) {
						d.set(i, temp * d.get(i) - g.get(i) - hs.get(i));
						final double t = d.get(i);
						dd += t * t;
						ds += t * step[i];
						ss += step[i] * step[i];
					}
					if (ds <= 0) {
						state = 160;
						continue;
					}
					if (ss < delSq) {
						state = 40;
						continue;
					}
				}
				crvMin[0] = 0;
				iterSw = iterC;
				
//				Test whether an alternative iteration is required.
				
			case 90:
				if (gg <= 1.0E-4 * ggBeg) {
					state = 160;
					continue;
				}
				sg = 0;
				shs = 0;
				for (int i = 0; i < n; i++) {
					sg += step[i] * g.get(i);
					shs += step[i] * hs.get(i);
				}
				final double sgk = sg + shs;
				final double angTest = sgk / Math.sqrt(gg * delSq);
				if (angTest <= -0.99) {
					state = 160;
					continue;
				}
				
//				Begin the alternative iteration by calculating D and HD and some
//				scalar products.

				iterC++;
				temp = Math.sqrt(delSq * gg - sgk * sgk);
				tempA = delSq / temp;
				tempB = sgk / temp;
				for (int i = 0; i < n; i++) {
					d.set(i, tempA * (g.get(i) + hs.get(i)) - tempB * step[i]);
				}
				state = 170;
				continue;
				
			case 120:
				double dg = 0;
				dhd = 0;
				double dhs = 0;
				for (int i = 0; i < n; i++) {
					dg += d.get(i) * g.get(i);
					dhd += hd.get(i) * d.get(i);
					dhs += hd.get(i) * step[i];
				}
							
//				Seek the value of the angle that minimizes Q.
				
				final double cf = 0.5 * (shs - dhd);
				final double qBeg = sg + cf;
				double qSav = qBeg;
				double qMin = qBeg;
				int iSave = -1;
				final int iu = 49;
				temp = TWO_PI / (double) (iu + 1);
				double qNew = 0;
				double angle;
				double cth;
				double sth;
				for (int i = 0; i < iu; i++) {
					angle = (double) (i + 1) * temp;
					cth = Math.cos(angle);
					sth = Math.sin(angle);
					qNew = (sg + cf * cth) * cth + (dg + dhs * cth) * sth;
					if (qNew < qMin) {
						qMin = qNew;
						iSave = i;
						tempA = qSav;
					} else if (i == iSave + 1) {
						tempB = qNew;
					}
					qSav = qNew;
				}
				if (iSave == -1.0) {
					tempA = qNew;
				}
				if (iSave == (iu - 1)) {
					tempB = qBeg;
				}
				angle = 0;
				if (tempA != tempB) {
					tempA -= qMin;
					tempB -= qMin;
					angle = 0.5 * (tempA - tempB) / (tempA + tempB);
				}
				angle = temp * ((double) (iSave + 1) + angle);
			
//				Calculate the new STEP and HS. Then test for convergence.
				
				cth = Math.cos(angle);
				sth = Math.sin(angle);
				final double reduc = qBeg - (sg + cf * cth) * cth - (dg + dhs * cth) * sth;
				gg = 0;
				for (int i = 0; i < n; i++) {
					step[i] = cth * step[i] + sth * d.get(i);
					hs.set(i, cth * hs.get(i) + sth * hd.get(i));
					final double t = g.get(i) + hs.get(i);
					gg += t * t;
				}
				qRed += reduc;
				final double ratio = reduc / qRed;
				if (iterC < iterMax && ratio >= 0.01) {
					state = 90;
					continue;
				}
				
			case 160:
				return;
								
//				The following instructions act as a subroutine for setting the vector
//				HD to the vector D multiplied by the second derivative matrix of Q.
//				They are called from three different places, which are distinguished
//				by the value of ITERC.
				
			case 170:
				for (int i = 0; i < n; i++) {
					hd.set(i, 0);
				}
				for (int k = 0; k < npt; k++) {
					temp = 0;
					for (int j = 0; j < n; j++) {
						temp += xPt[k][j] * d.get(j);
					}
					temp *= pq[k];
					for (int i = 0; i < n; i++) {
						hd.set(i, hd.get(i) + temp * xPt[k][i]);
					}
				}
				int ih = -1;
				for (int j = 0; j < n; j++) {
					for (int i = 0; i <= j; i++) {
						ih++;
						if (i < j) {
							hd.set(j, hd.get(j) + hq[ih] * d.get(i));
						}
						hd.set(i, hd.get(i) + hq[ih] * d.get(j));
					}
				}
				if (iterC == 0) {
					state = 20;
					continue;
				}
				if (iterC <= iterSw) {
					state = 50;
					continue;
				}
				state = 120;
				continue;				
			}
		}
	}
	
	
	/**
	  @param N is the number of variables.
      @param NPT is the number of interpolation equations.
      @param XOPT is the best interpolation point so far.
      @param XPT contains the coordinates of the current interpolation points.
      @param BMAT provides the last N columns of H.
      @param ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
      @param NDIM is the first dimension of BMAT and has the value NPT+N.
      @param KNEW is the index of the interpolation point that is going to be moved.
      @param DELTA is the current trust region bound.
      @param D will be set to the step from XOPT to the new point.
      @param ALPHA will be set to the KNEW-th diagonal element of the H matrix.
      @param HCOL, GC, GD, S and W will be used for working space.
 
      The step D is calculated in a way that attempts to maximize the modulus
      of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
      the KNEW-th Lagrange function
	 */
        @SuppressWarnings("fallthrough")
	private void biglag(final int n, final int npt, final double[] xOpt,
			final double[][] xpt, final double[][] bMat, final double[][] zMat,
			final int idz, final int nDim, final int knew, final double delta,
			final double[] d, final double[] alpha, final ScopedPtr hCol,
			final ScopedPtr gc, final ScopedPtr gd, final ScopedPtr s,
			final ScopedPtr w) {
		
//		Set some constants.

		final double delSq = delta * delta;
		final int nptm = npt - n - 1;
				
//		Set the first NPT components of HCOL to the leading elements of the
//		KNEW-th column of H.
		
		int iterC = -1;
		double temp;		
		for (int k = 0; k < npt; k++) {
			hCol.set(k, 0);			
		}
		for (int j = 0; j < nptm; j++) {
			temp = zMat[knew][j];
			if (j < idz - 1) {
				temp = -temp;
			}
			for (int k = 0; k < npt; k++) {
				hCol.set(k, hCol.get(k) + temp * zMat[k][j]);
			}			
		}
		alpha[0] = hCol.get(knew);

//		Set the unscaled initial direction D. Form the gradient of LFUNC at
//		XOPT, and multiply D by the second derivative matrix of LFUNC.

		double dd = 0;
		for (int i = 0; i < n; i++) {
			d[i] = xpt[knew][i] - xOpt[i];
			gc.set(i, bMat[knew][i]);
			gd.set(i, 0);
			dd += d[i] * d[i];
		}
		double sum;
		for (int k = 0; k < npt; k++) {
			temp = 0;
			sum = 0;
			for (int j = 0; j < n; j++) {
				temp += xpt[k][j] * xOpt[j];
				sum += xpt[k][j] * d[j];
			}
			temp *= hCol.get(k);
			sum *= hCol.get(k);
			for (int i = 0; i < n; i++) {
				gc.set(i, gc.get(i) + temp * xpt[k][i]);
				gd.set(i, gd.get(i) + sum * xpt[k][i]);
			}
		}
				
//		Scale D and GD, with a sign change if required. Set S to another
//		vector in the initial two dimensional subspace.

		double gg = 0;
		double sp = 0;
		double dhd = 0;
		for (int i = 0; i < n; i++) {
			final double t = gc.get(i);
			gg += t * t;
			sp += d[i] * t;
			dhd += d[i] * gd.get(i);
		}
		double scale = delta / Math.sqrt(dd);
		if (sp * dhd < 0) {
			scale = -scale;
		}
		temp = 0;
		if (sp * sp > 0.99 * dd * gg) {
			temp = 1;
		}
		double tau = scale * (Math.abs(sp) + 0.5 * scale * Math.abs(dhd));
		if (gg * delSq < 0.01 * tau * tau) {
			temp = 1;
		}
		for (int i = 0; i < n; i++) {
			d[i] *= scale;
			gd.set(i, scale * gd.get(i));
			s.set(i, gc.get(i) + temp * gd.get(i));
		}
		
//		Begin the iteration by overwriting S with a vector that has the
//		required length and direction, except that termination occurs if
//		the given D and S are nearly parallel.
		
		int state = 80;
		for (;;) {
			switch (state) {
			case 80:
				iterC++;
				dd = 0;
				sp = 0;
				double ss = 0;
				for (int i = 0; i < n; i++) {
					dd += d[i] * d[i];
					final double t = s.get(i);
					sp += d[i] * t;
					ss += t * t;
				}
				temp = dd * ss - sp * sp;
				if (temp <= 1E-8 * dd * ss) {
					state = 160;
					continue;
				}
				double denom = Math.sqrt(temp);
				for (int i = 0; i < n; i++) {
					s.set(i, (dd * s.get(i) - sp * d[i]) / denom);
					w.set(i, 0);
				}
								
//				Calculate the coefficients of the objective function on the circle,
//				beginning with the multiplication of S by the second derivative matrix.
				
				for (int k = 0; k < npt; k++) {
					sum = 0;
					for (int j = 0; j < n; j++) {
						sum += xpt[k][j] * s.get(j);
					}
					sum *= hCol.get(k);
					for (int i = 0; i < n; i++) {
						w.set(i, w.get(i) + sum * xpt[k][i]);
					}
				}
				double cf1 = 0;
				double cf2 = 0;
				double cf3 = 0;
				double cf4 = 0;
				double cf5 = 0;
				for (int i = 0; i < n; i++) {
					cf1 += s.get(i) * w.get(i);
					cf2 += d[i] * gc.get(i);
					cf3 += s.get(i) * gc.get(i);
					cf4 += d[i] * gd.get(i);
					cf5 += s.get(i) * gd.get(i);
				}
				cf1 *= 0.5;
				cf4 = 0.5 * cf4 - cf1;

//				Seek the value of the angle that maximizes the modulus of TAU.

				double tauBeg = cf1 + cf2 + cf4;
				double tauMax = tauBeg;
				double tauOld = tauBeg;
				double tempA = 0;
				double tempB = 0;
				double angle = 0;
				double cth;
				double sth;
				int iSave = -1;				
				int iu = 49;
				temp = TWO_PI / (double) (iu + 1);
				for (int i = 0; i < iu; i++) {
					angle = (double) (i + 1) * temp;
					cth = Math.cos(angle);
					sth = Math.sin(angle);
					tau = cf1 + (cf2 + cf4 * cth) * cth + (cf3 + cf5 * cth) * sth;
					if (Math.abs(tau) > Math.abs(tauMax)) {
						tauMax = tau;
						iSave = i;
						tempA = tauOld;
					} else if (i == iSave + 1) {
						tempB = tau;
					}
					tauOld = tau;
				}
				if (iSave == -1) {
					tempA = tau;
				}
				if (iSave == (iu - 1)) {
					tempB = tauBeg;
				}
				double step = 0;
				if (tempA != tempB) {
					tempA -= tauMax;
					tempB -= tauMax;
					step = 0.5 * (tempA - tempB) / (tempA + tempB);
				}
				angle = temp * ((double) (iSave + 1) + step);

//				Calculate the new D and GD. Then test for convergence.
				
				cth = Math.cos(angle);
				sth = Math.sin(angle);
				tau = cf1 + (cf2 + cf4 * cth) * cth + (cf3 + cf5 * cth)* sth;
				for (int i = 0; i < n; i++) {
					d[i] = cth * d[i] + sth * s.get(i);
					gd.set(i, cth * gd.get(i) + sth * w.get(i));
					s.set(i, gc.get(i) + gd.get(i));
				}
				if (Math.abs(tau) <= 1.1 * Math.abs(tauBeg)) {
					state = 160;
					continue;
				}
				if (iterC < n) {
					state = 80;
					continue;
				}
				
			case 160:
				return;
			}
		}
	}
	
	
	
	/**
	 @param N is the number of variables.
     @param NPT is the number of interpolation equations.
     @param XOPT is the best interpolation point so far.
     @param XPT contains the coordinates of the current interpolation points.
     @param BMAT provides the last N columns of H.
     @param ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
     @param NDIM is the first dimension of BMAT and has the value NPT+N.
     @param KOPT is the index of the optimal interpolation point.
     @param KNEW is the index of the interpolation point that is going to be moved.
     @param D will be set to the step from XOPT to the new point, and on entry it
       should be the D that was calculated by the last call of BIGLAG. The
       length of the initial D provides a trust region bound on the final D.
     @param W will be set to Wcheck for the final choice of D.
     @param VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
     @param BETA will be set to the value that will occur in the updating formula
       when the KNEW-th interpolation point is moved to its new position.
     @param S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
       for working space.

     D is calculated in a way that should provide a denominator with a large
     modulus in the updating formula when the KNEW-th interpolation point is
     shifted to the new position XOPT+D.
	 */
        @SuppressWarnings("fallthrough")
	private void bigDen(final int n, final int npt, final double[] xOpt,
			final double[][] xpt, final double[][] bMat, final double[][] zMat,
			final int idz, final int nDim, final int kOpt, final int kNew,
			final double[] d, final ScopedPtr w, final double[] vLag, double beta,
			final double[] s, final double[][] wVec, final double[][] prod) {
		
//	    		      DIMENSION DEN(9),DENEX(9),PAR(9)
		final double[] den = new double[9];
		final double[] denex = new double[9];
		final double[] par = new double[9];
		
//	    		C     Set some constants.

		final int nptm = npt - n - 1;

//		Store the first NPT elements of the KNEW-th column of H in W(N+1)
//		to W(N+NPT).
		
		for (int k = 0; k < npt; k++) {
			w.set(n + k, 0);
		}
		double temp;
		for (int j = 0; j < nptm; j++) {
			temp = zMat[kNew][j];
			if (j < idz - 1) {
				temp = -temp;
			}
			for (int k = 0; k < npt; k++) {
				w.set(n + k, w.get(n + k) + temp * zMat[k][j]);
			}
		}
		double alpha = w.get(n + kNew);
		
//		The initial search direction D is taken from the last call of BIGLAG,
//		and the initial S is set below, usually to the direction from X_OPT
//		to X_KNEW, but a different direction to an interpolation point may
//		be chosen, in order to prevent S from being nearly parallel to D.

		double dd = 0;
		double ds = 0;
		double ss = 0;
		double xOptSq = 0;
		for (int i = 0; i < n; i++) {
			dd += d[i] * d[i];
			s[i] = xpt[kNew][i] - xOpt[i];
			ds += d[i] * s[i];
			ss += s[i] * s[i];
			xOptSq += xOpt[i] * xOpt[i];
		}
		if (ds * ds > 0.99 * dd * ss) {
			int kSav = kNew;
			double dTest = ds * ds / ss;
			for (int k = 0; k < npt; k++) {
				if (k != kOpt) {
					double dsTemp = 0;
					double ssTemp = 0;
					for (int i = 0; i < n; i++) {
						double diff = xpt[k][i] - xOpt[i];
						dsTemp += d[i] * diff;
						ssTemp += diff * diff;
					}
					if (dsTemp * dsTemp / ssTemp < dTest) {
						kSav = k;
						dTest = dsTemp * dsTemp / ssTemp;
						ds = dsTemp;
						ss = ssTemp;
					}
				}
			}
			for (int i = 0; i < n; i++) {
				s[i] = xpt[kSav][i] - xOpt[i];			
			}
		}
		double ssDen = dd * ss - ds * ds;
		int iterC = 0;
		double denSav = 0;
		
//		Begin the iteration by overwriting S with a vector that has the
//		required length and direction.
		
		double xOptD = 0;
		double xOptS = 0;
		double sum;
		double tempC;
		
		int state = 70;
		for (;;) {
			switch (state) {
			case 70:
				iterC++;
				temp = 1.0 / Math.sqrt(ssDen);
				xOptD = 0;
				xOptS = 0;
				for (int i = 0; i < n; i++) {
					s[i] = temp * (dd * s[i] - ds * d[i]);
					xOptD += xOpt[i] * d[i];
					xOptS += xOpt[i] * s[i];
				}
				
//				Set the coefficients of the first two terms of BETA.				
				
				double tempA = 0.5 * xOptD * xOptD;
				double tempB = 0.5 * xOptS * xOptS;
				den[0] = dd * (xOptSq + 0.5 * dd) + tempA + tempB;
				den[1] = 2.0 * xOptD * dd;
				den[2] = 2.0 * xOptS * dd;
				den[3] = tempA - tempB;
				den[4] = xOptD * xOptS;
				for (int i = 5; i < 9; i++) {
					den[i] = 0;
				}

//				Put the coefficients of Wcheck in WVEC.				
				
				for (int k = 0; k < npt; k++) {
					tempA = 0;
					tempB = 0;
					tempC = 0;
					for (int i = 0; i < n; i++) {
						tempA += xpt[k][i] * d[i];
						tempB += xpt[k][i] * s[i];
						tempC += xpt[k][i] * xOpt[i];
					}
					wVec[k][1] = 0.25 * (tempA * tempA + tempB * tempB);
					wVec[k][2] = tempA * tempC; 
					wVec[k][3] = tempB * tempC;
					wVec[k][4] = 0.25 * (tempA * tempA - tempB * tempB);
					wVec[k][5] = 0.5 * tempA * tempB;
				}
				for (int i = 0; i < n; i++) {
					int ip = i + npt;
					wVec[ip][1] = 0;
					wVec[ip][2] = d[i];
					wVec[ip][3] = s[i];
					wVec[ip][4] = 0;
					wVec[ip][5] = 0;
				}

//				Put the coefficents of THETA*Wcheck in PROD.
				
				for (int jc = 0; jc < 5; jc++) {
					int nw = npt;
					if (jc == 1 || jc == 2) {
						nw = nDim;
					}
					for (int k = 0; k < npt; k++) {
						prod[k][jc] = 0;
					}
					for (int j = 0; j < nptm; j++) {
						sum = 0;
						for (int k = 0; k < npt; k++) {
							sum += zMat[k][j] * wVec[k][jc];
						}
						if (j < idz - 1) {
							sum = -sum;
						}
						for (int k = 0; k < npt; k++) {
							prod[k][jc] += sum * zMat[k][j];
						}
					}
					if (nw == nDim) {
						for (int k = 0; k < npt; k++) {
							sum = 0;
							for (int j = 0; j < n; j++) {
								sum += bMat[k][j] * wVec[npt + j][jc];
								prod[k][jc] += sum;
							}
						}
					}
					for (int j = 0; j < n; j++) {
						sum = 0;
						for (int i = 0; i < nw; i++) {
							sum += bMat[i][j] * wVec[i][jc];
						}
						prod[npt + j][jc] = sum;
					}
				}
								
//				Include in DEN the part of BETA that depends on THETA.
				
				for (int k = 0; k < nDim; k++) {
					sum = 0;
					for (int i = 0; i < 5; i++) {
						par[i] = 0.5 * prod[k][i] * wVec[k][i];
						sum += par[i];
					}
					den[0] -= par[0] - sum;
					tempA = prod[k][0] * wVec[k][1] + prod[k][1] * wVec[k][0];
					tempB = prod[k][1] * wVec[k][3] + prod[k][3] * wVec[k][1];
					tempC = prod[k][2] * wVec[k][4] + prod[k][4] * wVec[k][2];
					den[1] -= tempA - 0.5 * (tempB + tempC);
					den[5] -= 0.5 * (tempB - tempC);
					tempA = prod[k][0] * wVec[k][2] + prod[k][2] * wVec[k][0];
					tempB = prod[k][1] * wVec[k][4] + prod[k][4] * wVec[k][1];
					tempC = prod[k][2] * wVec[k][3] + prod[k][3] * wVec[k][2];
					den[2] -= tempA - 0.5 * (tempB - tempC);
					den[6] -= 0.5 * (tempB + tempC);
					tempA = prod[k][0] * wVec[k][3] + prod[k][3] * wVec[k][0];
					den[3] -= tempA - par[1] + par[2];
					tempA = prod[k][0] * wVec[k][4] + prod[k][4] * wVec[k][0];
					tempB = prod[k][1] * wVec[k][2] + prod[k][2] * wVec[k][1];
					den[4] -= tempA - 0.5 * tempB;
					den[7] -= par[3] + par[4];
					tempA = prod[k][3] * wVec[k][4] + prod[k][4] * wVec[k][3];
					den[8] -= 0.5 * tempA;
				}

//				Extend DEN so that it holds all the coefficients of DENOM.

				sum = 0;
				for (int i = 0; i < 5; i++) {
					par[i] = 0.5 * prod[kNew][i] * prod[kNew][i];
					sum += par[i];
				}
				denex[0] = alpha * den[0] + par[0] + sum;
				tempA = 2.0 * prod[kNew][0] * prod[kNew][1];
				tempB = prod[kNew][1] * prod[kNew][3];
				tempC = prod[kNew][2] * prod[kNew][4];
				denex[1] = alpha * den[1] + tempA + tempB + tempC;
				denex[5] = alpha * den[5] + tempB - tempC;
				tempA = 2.0 * prod[kNew][0] * prod[kNew][2];
				tempB = prod[kNew][1] * prod[kNew][4];
				tempC = prod[kNew][2] * prod[kNew][3];
				denex[2] = alpha * den[2] + tempA + tempB - tempC;
				denex[6] = alpha * den[6] + tempB + tempC;
				tempA = 2.0 * prod[kNew][0] * prod[kNew][3];
				denex[3] = alpha * den[3] + tempA + par[1] - par[2];
				tempA = 2.0 * prod[kNew][0] * prod[kNew][4];
				denex[4] = alpha * den[4] + tempA + prod[kNew][1] * prod[kNew][2];
				denex[7] = alpha * den[7] + par[3] - par[4];
				denex[8] = alpha * den[8] + prod[kNew][3] * prod[kNew][4];

//				Seek the value of the angle that maximizes the modulus of DENOM.

				sum = denex[0] + denex[1] + denex[3] + denex[5] + denex[7];
				double denOld = sum;
				double denMax = sum;
				int iSave = 0;
				int iu = 49;
				temp = TWO_PI / (double) (iu + 1);
				par[0] = 1;
				double angle;
				for (int i = 0; i < iu; i++) {
					angle = (double) i * temp;
					par[1] = Math.cos(angle);
					par[2] = Math.sin(angle);
					for (int j = 3; j < 8; j = j + 2) {
						par[j] = par[1] * par[j - 2] - par[2] * par[j - 1];
						par[j + 1] = par[1] * par[j - 1] + par[2] * par[j - 2];
					}
					double sumOld = sum;
					sum = 0;
					for (int j = 0; j < 9; j++) {
						sum += denex[j] * par[j];
					}
					if (Math.abs(sum) > Math.abs(denMax)) {
						denMax = sum;
						iSave = i;
						tempA = sumOld;
					} else if (i == iSave + 1) {
						tempB = sum;
					}
				}
				if (iSave == 0) {
					tempA = sum;
				}
				if (iSave == iu) {
					tempB = denOld;
				}
				double step = 0;
				if (tempA != tempB) {
					tempA -= denMax;
					tempB -= denMax;
					step = 0.5 * (tempA - tempB) / (tempA + tempB);
				}
				angle = temp * ((double) iSave + step);
				
//				 Calculate the new parameters of the denominator, the new VLAG vector
//				 and the new D. Then test for convergence.
				
				par[1] = Math.cos(angle);
				par[2] = Math.sin(angle);
				for (int j = 3; j < 8; j = j + 2) {
					par[j] = par[1] * par[j - 2] - par[2] * par[j - 1];
					par[j + 1] = par[1] * par[j - 1] + par[2] * par[j - 2];
				}
				beta = 0;
				denMax = 0;
				for (int j = 0; j < 9; j++) {
					beta += den[j] * par[j];
					denMax += denex[j] * par[j];
				}
				for (int k = 0; k < nDim; k++) {
					vLag[k] = 0;
					for (int j = 0; j < 5; j++) {
						vLag[k] += prod[k][j] * par[j];
					}
				}
				double tau = vLag[kNew];
				dd = 0;
				tempA = 0;
				tempB = 0;
				for (int i = 0; i < n; i++) {
					d[i] *= par[1] + par[2] * s[i];
					w.set(i, xOpt[i] + d[i]);
					dd += d[i] * d[i];
					final double t = w.get(i);
					tempA += d[i] * t;
					tempB += t * t;
				}
				if (iterC >= n) {
					state = 340;
					continue;
				}
				if (iterC > 0) {
					denSav = Math.max(denSav, denOld);
				}
				if (Math.abs(denMax) <= 1.1 * Math.abs(denSav)) {
					state = 340;
					continue;
				}
				denSav = denMax;
				
//				Set S to half the gradient of the denominator with respect to D.
//				Then branch for the next iteration.
				
				for (int i = 0; i < n; i++) {
					temp = tempA * xOpt[i] + tempB * d[i] - vLag[npt + i];
					s[i] = tau * bMat[kNew][i] + alpha * temp;
				}
				for (int k = 0; k < npt; k++) {
					sum = 0;
					for (int j = 0; j < n; j++) {
						sum += xpt[k][j] * w.get(j);
					}
					temp = (tau * w.get(n + k) - alpha * vLag[k]) * sum;
					for (int i = 0; i < n; i++) {
						s[i] += temp * xpt[k][i];
					}
				}
				ss = 0;
				ds = 0;
				for (int i = 0; i < n; i++) {
					ss += s[i] * s[i];
					ds += d[i] * s[i];
				}
				ssDen = dd * ss - ds * ds;
				if (ssDen >= 1E-8 * dd * ss) {
					state = 70;
					continue;
				}
				
//				Set the vector W before the RETURN from the subroutine.
				
			case 340:
				for (int k = 0; k < nDim; k++) {
					w.set(k, 0);
					for (int j = 0; j < 5; j++) {
						w.set(k, w.get(k) + wVec[k][j] * par[j]);
					}
				}
				vLag[kOpt] += 1;
				return;
			}
		}
	}
	
	/**
	 * The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
      interpolation point that has index KNEW. On entry, VLAG contains the
      components of the vector Theta*Wcheck+e_b of the updating formula
      (6.11), and BETA holds the value of the parameter that has this name.
      The vector W is used for working space.
	 */
	private void update(final int n, final int npt, final double[][] bMat, final double[][] zMat,
			int idz, final int nDim, final double[] vLag, final double beta, final int kNew,
			final ScopedPtr w) {
		
		final int nptm = npt - n - 1;

//		     Apply the rotations that put zeros in the KNEW-th row of ZMAT.

		double temp;
		double tempA;
		double tempB = 0;
		int jl = 0;
		for (int j = 1; j < nptm; j++) {
			if (j == idz - 1) {
				jl = idz - 1;
			} else if (zMat[kNew][j] != 0) {
				temp = Math.sqrt(zMat[kNew][jl] * zMat[kNew][jl] + zMat[kNew][j] * zMat[kNew][j]);
				tempA = zMat[kNew][jl] / temp;
				tempB = zMat[kNew][j] / temp;
				for (int i = 0; i < npt; i++) {
					temp = tempA * zMat[i][jl] + tempB * zMat[i][j];
					zMat[i][j] = tempA * zMat[i][j] - tempB * zMat[i][jl];
					zMat[i][jl] = temp;
				}
				zMat[kNew][j] = 0;
			}
		}
		
//		Put the first NPT components of the KNEW-th column of HLAG into W,
//		and calculate the parameters of the updating formula.
		
		tempA = zMat[kNew][0];
		if (idz > 1) { 
			tempA = -tempA;
		}
		if (jl > 0) {
			tempB = zMat[kNew][jl];
		}
		for (int i = 0; i < npt; i++) {
			w.set(i, tempA * zMat[i][0]);
			if (jl > 0) {
				w.set(i, w.get(i) + tempB * zMat[i][jl]);
			}
		}		
		double alpha = w.get(kNew);
		double tau = vLag[kNew];
		double tauSq = tau * tau;
		double denom = alpha * beta + tauSq;
		vLag[kNew] -= 1;
		
//		Complete the updating of ZMAT when there is only one nonzero element
//		in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
//		then the first column of ZMAT will be exchanged with another one later.
		
		int iFlag = 0;
		if (jl == 0) {
			temp = Math.sqrt(Math.abs(denom));
			tempB = tempA / temp;
			tempA = tau / temp;
			for (int i = 0; i < npt; i++) {
				zMat[i][0] = zMat[i][0] * tempA - tempB * vLag[i];
			}
			if (idz == 1 && temp < 0) {
				idz = 2;
			}
			if (idz >= 2 && temp >= 0) {
				iFlag = 1;
			}
		} else {
			
//			Complete the updating of ZMAT in the alternative case.
			
			int ja = 0;
			if (beta >= 0) {
				ja = jl;
			}
			int jb = jl + 1 - ja;
			temp = zMat[kNew][jb] / denom;
			tempA = temp * beta;
			tempB = temp * tau;
			temp = zMat[kNew][ja];
			double scalA = 1.0 / Math.sqrt(Math.abs(beta) * temp * temp + tauSq);
			double scalB = scalA * Math.sqrt(Math.abs(denom));
			for (int i = 0; i < npt; i++) {
				zMat[i][ja] = scalA * (tau * zMat[i][ja] - temp * vLag[i]);
				zMat[i][jb] = scalB * (zMat[i][jb] - tempA * w.get(i) - tempB * vLag[i]);
			}
			if (denom <= 0) {
				if (beta < 0) {
					idz++;
				}
				if (beta >= 0) {
					iFlag = 1;
				}
			}
		}
		
//		IDZ is reduced in the following case, and usually the first column
//		of ZMAT is exchanged with a later one.
		
		if (iFlag == 1) {
			idz--;
			for (int i = 0; i < npt; i++) {
				temp = zMat[i][0];
				zMat[i][0] = zMat[i][idz];
				zMat[i][idz] = temp;
			}
		}
		
//		Finally, update the matrix BMAT.
		
		for (int j = 0; j < n; j++) {
			int jp = npt + j;
			w.set(jp, bMat[kNew][j]);
			tempA = (alpha * vLag[jp] - tau * w.get(jp)) / denom;
			tempB = (-beta * w.get(jp) - tau * vLag[jp]) / denom;
			for (int i = 0; i <= jp; i++) {
				bMat[i][j] += tempA * vLag[i] + tempB * w.get(i);
				if (i >= npt) {
					bMat[jp][i-npt] = bMat[i][j];
				}
			}
		}
		return;
	}
	
    /**
     * Used to simulate Fortran pointers.
     * 
     * XXX optimization: implement +=
     */
    private static class ScopedPtr {
        /**
         * array storing elements.
         */
        private double[] w;
        /**
         * base index for access.
         */
        private int base;

        /**
         * @param w array storing elements.
         * @param base base index for access.
         */
        ScopedPtr(double[] w, int base) {
            this.w = w;
            this.base = base;
        }

        /**
         * @param index relative index of returned ScopedPtr
         * @return ScopedPtr with new base = this.base + index
         */
        ScopedPtr ptr(int index) {
            return new ScopedPtr(w, base + index);
        }

        /**
         * @param index of accessed element relative to base.
         * @return value returned value at index.
         */
        double get(int index) {
            return w[base + index];
        }

//        /**
//         * @param index of accessed elements relative to base.
//         * @param n number of values to be returned.
//         * @return n values starting at index.
//         */
//        double[] getAll(int index, int n) {
//            return Arrays.copyOfRange(w, base+index, base+index+n);
//        }

//        /**
//         * @return all elements.
//         */
//		double[] getAll() {
//            return w;
//        }

        /**
         * @param index index of accessed element relative to base.
         * @param value stored at index.
         */
        void set(int index, double value) {
            w[base + index] = value;
        }

        /* (non-Javadoc)
         * @see java.lang.Object#toString()
         */
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < 20; i++)
                if (base + i >= 0 && base + i < w.length)
                    sb.append("").append(i).append(":").append(w[base + i]).append("\n");
            return sb.toString();
        }

    }
}
