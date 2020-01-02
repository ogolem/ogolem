/* RISO: an implementation of distributed belief networks.
 * Copyright (C) 1999, Robert Dodier.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA,
 * or visit the GNU web site, www.gnu.org.
 */
package contrib.lbfgs;

import java.io.Serializable;

/** <p> This class contains code for the limited-memory Broyden-Fletcher-Goldfarb-Shanno
  * (LBFGS) algorithm for large-scale multidimensional unconstrained minimization problems.
  * This file is a translation of Fortran code written by Jorge Nocedal.
  * The only modification to the algorithm is the addition of a cache to 
  * store the result of the most recent line search. See <var>solution_cache</var> below.
  *
  * LBFGS is distributed as part of the RISO project. Following is a message from Jorge Nocedal:
  * <pre>
  *   From: Jorge Nocedal [mailto:nocedal@dario.ece.nwu.edu]
  *   Sent: Friday, August 17, 2001 9:09 AM
  *   To: Robert Dodier
  *   Subject: Re: Commercial licensing terms for LBFGS?
  *   
  *   Robert:
  *   The code L-BFGS (for unconstrained problems) is in the public domain.
  *   It can be used in any commercial application.
  *   
  *   The code L-BFGS-B (for bound constrained problems) belongs to
  *   ACM. You need to contact them for a commercial license. It is
  *   algorithm 778.
  *   
  *   Jorge
  * </pre>
  * 
  * <p> This code is derived from the Fortran program <code>lbfgs.f</code>.
  * The Java translation was effected mostly mechanically, with some
  * manual clean-up; in particular, array indices start at 0 instead of 1.
  * Most of the comments from the Fortran code have been pasted in here
  * as well.</p>
  *
  * <p> Here's some information on the original LBFGS Fortran source code,
  * available at <a href="http://www.netlib.org/opt/lbfgs_um.shar">
  * http://www.netlib.org/opt/lbfgs_um.shar</a>. This info is taken
  * verbatim from the Netlib blurb on the Fortran source.</p>
  *
  * <pre>
  * 	file    opt/lbfgs_um.shar
  * 	for     unconstrained optimization problems
  * 	alg     limited memory BFGS method
  * 	by      J. Nocedal
  * 	contact nocedal@eecs.nwu.edu
  * 	ref     D. C. Liu and J. Nocedal, ``On the limited memory BFGS method for
  * 	,       large scale optimization methods'' Mathematical Programming 45
  * 	,       (1989), pp. 503-528.
  * 	,       (Postscript file of this paper is available via anonymous ftp
  * 	,       to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_um.)
  * </pre>
  *
  * @author Jorge Nocedal: original Fortran version, including comments
  * (July 1990). Robert Dodier: Java translation, August 1997.
  * @author Johannes Dieterich: some performance optimizations in relevant parts
  * (June 2012)
  * @author Johannes Dieterich: some stabilization patches in the linesearch return
  * Also trying to make the output to console less confusing.
  * (April 2013)
  * @author Johannes Dieterich: marking the inner class static and fixing some
  * errors that must have been oversights by the original author (stpmin to  stepmin)
  * and are caused by having everything in class variables.
  * (January 2014)
  */

public class LBFGS implements Serializable {
    
    private static final long serialVersionUID = (long) 20111117;
        
    /** Specialized exception class for LBFGS; contains the
     * <code>iflag</code> value returned by <code>lbfgs</code>.
     */
    public static class ExceptionWithIflag extends Exception {

        private static final long serialVersionUID = 1L;
        public int iflag;

        public ExceptionWithIflag(int i, String s) {
            super(s);
            iflag = i;
        }

        @Override
        public String toString() {
            return getMessage() + " (iflag == " + iflag + ")";
        }
    }
    /** Controls the accuracy of the line search <code>mcsrch</code>. If the
     * function and gradient evaluations are inexpensive with respect
     * to the cost of the iteration (which is sometimes the case when
     * solving very large problems) it may be advantageous to set <code>gtol</code>
     * to a small value. A typical small value is 0.1.  Restriction:
     * <code>gtol</code> should be greater than 1e-4.
     */
    public final double gtol;
    /** Specify lower bound for the step in the line search.
     * The default value is 1e-20. This value need not be modified unless
     * the exponent is too large for the machine being used, or unless
     * the problem is extremely badly scaled (in which case the exponent
     * should be increased).
     */
    public final double stpmin = 1e-20;
    /** Specify upper bound for the step in the line search.
     * The default value is 1e20. This value need not be modified unless
     * the exponent is too large for the machine being used, or unless
     * the problem is extremely badly scaled (in which case the exponent
     * should be increased).
     */
    public final double stpmax = 1e20;
    /** The solution vector as it was at the end of the most recently
     * completed line search. This will usually be different from the
     * return value of the parameter <var>x</var> of <var>lbfgs</var>, which
     * is modified by line-search steps. A caller which wants to stop the
     * optimization iterations before <var>LBFGS.lbfgs</var> automatically stops
     * (by reaching a very small gradient) should copy this vector instead
     * of using <var>x</var>. When <var>LBFGS.lbfgs</var> automatically stops,
     * then <var>x</var> and <var>solution_cache</var> are the same.
     */
    public double[] solution_cache = null;
    private double gnorm = 0, stp1 = 0, ftol = 0, stp[] = new double[1], ys = 0, yy = 0, sq = 0, yr = 0, beta = 0, xnorm = 0;
    private int iter = 0, nfun = 0, point = 0, ispt = 0, iypt = 0, info[] = new int[1], bound = 0, npt = 0, cp = 0, nfev[] = new int[1], inmc = 0, iycn = 0, iscn = 0;
    private boolean finish = false;
    private double[] w = null;

    // Mcsrch (is a linesearch)
    private final Mcsrch mcsrch;
    private final int maxfev;
    
    public LBFGS(final double linetol, final int maxIterLine) {
        if (linetol <= 0.0001) {
            System.err.println("WARNING: LBFGS gtol is less than or equal to 0.0001. It has been reset to 0.9.");
            gtol = 0.9;
        } else gtol = linetol;

        this.maxfev = maxIterLine;
        this.mcsrch = new Mcsrch(gtol, stpmin, stpmax);
    }
    
    public LBFGS(final LBFGS orig){
        this.gtol = orig.gtol;
        this.maxfev = orig.maxfev;
        this.mcsrch = new Mcsrch(orig.mcsrch);
    }

    /** This method returns the total number of evaluations of the objective
     * function since the last time LBFGS was restarted. The total number of function
     * evaluations increases by the number of evaluations required for the
     * line search; the total is only increased after a successful line search.
     * @return number of function evaluations
     */
    public int nfevaluations() {
        return nfun;
    }

    /** This subroutine solves the unconstrained minimization problem
     * <pre>
     *     min f(x),    x = (x1,x2,...,x_n),
     * </pre>
     * using the limited-memory BFGS method. The routine is especially
     * effective on problems involving a large number of variables. In
     * a typical iteration of this method an approximation <code>Hk</code> to the
     * inverse of the Hessian is obtained by applying <code>m</code> BFGS updates to
     * a diagonal matrix <code>Hk0</code>, using information from the previous M steps.
     * The user specifies the number <code>m</code>, which determines the amount of
     * storage required by the routine. The user may also provide the
     * diagonal matrices <code>Hk0</code> if not satisfied with the default choice.
     * The algorithm is described in "On the limited memory BFGS method
     * for large scale optimization", by D. Liu and J. Nocedal,
     * Mathematical Programming B 45 (1989) 503-528.
     *
     * The user is required to calculate the function value <code>f</code> and its
     * gradient <code>g</code>. In order to allow the user complete control over
     * these computations, reverse  communication is used. The routine
     * must be called repeatedly under the control of the parameter
     * <code>iflag</code>. 
     *
     * The steplength is determined at each iteration by means of the
     * line search routine <code>mcsrch</code>, which is a slight modification of
     * the routine <code>CSRCH</code> written by More' and Thuente.
     *
     * The only variables that are machine-dependent are <code>xtol</code>,
     * <code>stpmin</code> and <code>stpmax</code>.
     *
     * Progress messages and non-fatal error messages are printed to <code>System.err</code>.
     * Fatal errors cause exception to be thrown, as listed below.
     *
     * @param n The number of variables in the minimization problem.
     *		Restriction: <code>n &gt; 0</code>.
     *
     * @param m The number of corrections used in the BFGS update. 
     *		Values of <code>m</code> less than 3 are not recommended;
     *		large values of <code>m</code> will result in excessive
     *		computing time. <code>3 &lt;= m &lt;= 7</code> is recommended.
     *		Restriction: <code>m &gt; 0</code>.
     *
     * @param x On initial entry this must be set by the user to the values
     *		of the initial estimate of the solution vector. On exit with
     *		<code>iflag = 0</code>, it contains the values of the variables
     *		at the best point found (usually a solution).
     *
     * @param f Before initial entry and on a re-entry with <code>iflag = 1</code>,
     *		it must be set by the user to contain the value of the function
     *		<code>f</code> at the point <code>x</code>.
     *
     * @param g Before initial entry and on a re-entry with <code>iflag = 1</code>,
     *		it must be set by the user to contain the components of the
     *		gradient <code>g</code> at the point <code>x</code>.
     *
     * @param diagco  Set this to <code>true</code> if the user  wishes to
     *		provide the diagonal matrix <code>Hk0</code> at each iteration.
     *		Otherwise it should be set to <code>false</code> in which case
     *		<code>lbfgs</code> will use a default value described below. If
     *		<code>diagco</code> is set to <code>true</code> the routine will
     *		return at each iteration of the algorithm with <code>iflag = 2</code>,
     *		and the diagonal matrix <code>Hk0</code> must be provided in
     *		the array <code>diag</code>.
     *
     * @param diag If <code>diagco = true</code>, then on initial entry or on
     *		re-entry with <code>iflag = 2</code>, <code>diag</code>
     *		must be set by the user to contain the values of the 
     *		diagonal matrix <code>Hk0</code>. Restriction: all elements of
     *		<code>diag</code> must be positive.
     *
     * @param iprint Specifies output generated by <code>lbfgs</code>.
     *		<code>iprint[0]</code> specifies the frequency of the output:
     *		<ul>
     *		<li> <code>iprint[0] &lt; 0</code>: no output is generated,
     *		<li> <code>iprint[0] = 0</code>: output only at first and last iteration,
     *		<li> <code>iprint[0] &gt; 0</code>: output every <code>iprint[0]</code> iterations.
     *		</ul>
     *
     *		<code>iprint[1]</code> specifies the type of output generated:
     *		<ul>
     *		<li> <code>iprint[1] = 0</code>: iteration count, number of function 
     *			evaluations, function value, norm of the gradient, and steplength,
     *		<li> <code>iprint[1] = 1</code>: same as <code>iprint[1]=0</code>, plus vector of
     *			variables and  gradient vector at the initial point,
     *		<li> <code>iprint[1] = 2</code>: same as <code>iprint[1]=1</code>, plus vector of
     *			variables,
     *		<li> <code>iprint[1] = 3</code>: same as <code>iprint[1]=2</code>, plus gradient vector.
     *		</ul>
     *
     *	@param eps Determines the accuracy with which the solution
     *		is to be found. The subroutine terminates when
     *		<pre>
     *            ||G|| &lt; EPS max(1,||X||),
     *		</pre>
     *		where <code>||.||</code> denotes the Euclidean norm.
     *
     *	@param xtol An estimate of the machine precision (e.g. 10e-16 on a
     *		SUN station 3/60). The line search routine will terminate if the
     *		relative width of the interval of uncertainty is less than
     *		<code>xtol</code>.
     *
     * @param iflag This must be set to 0 on initial entry to <code>lbfgs</code>.
     *		A return with <code>iflag &lt; 0</code> indicates an error,
     *		and <code>iflag = 0</code> indicates that the routine has
     *		terminated without detecting errors. On a return with
     *		<code>iflag = 1</code>, the user must evaluate the function
     *		<code>f</code> and gradient <code>g</code>. On a return with
     *		<code>iflag = 2</code>, the user must provide the diagonal matrix
     *		<code>Hk0</code>.
     *
     *		The following negative values of <code>iflag</code>, detecting an error,
     *		are possible:
     *		<ul>
     *		<li> <code>iflag = -1</code> The line search routine
     *			<code>mcsrch</code> failed. One of the following messages
     *			is printed:
     *			<ul>
     *			<li> Improper input parameters.
     *			<li> Relative width of the interval of uncertainty is at
     *				most <code>xtol</code>.
     *			<li> More than XX function evaluations were required at the
     *				present iteration.
     *			<li> The step is too small.
     *			<li> The step is too large.
     *			<li> Rounding errors prevent further progress. There may not
     *				be  a step which satisfies the sufficient decrease and
     *				curvature conditions. Tolerances may be too small.
     *			</ul>
     *		<li><code>iflag = -2</code> The i-th diagonal element of the diagonal inverse
     *			Hessian approximation, given in DIAG, is not positive.
     *		<li><code>iflag = -3</code> Improper input parameters for LBFGS
     *			(<code>n</code> or <code>m</code> are not positive).
     *		</ul>
     *
     *	@throws LBFGS.ExceptionWithIflag 
     */
    public void lbfgs(int n, int m, double[] x, double f, double[] g, boolean diagco, double[] diag, int[] iprint, double eps, double xtol, int[] iflag, boolean beReallySilent) throws ExceptionWithIflag {
        boolean execute_entire_while_loop = false;

        if (w == null || w.length != n * (2 * m + 1) + 2 * m) {
            w = new double[n * (2 * m + 1) + 2 * m];
        }

        if (iflag[0] == 0) {
            // Initialize.

            solution_cache = new double[n];
            System.arraycopy(x, 0, solution_cache, 0, n);

            iter = 0;

            if (n <= 0 || m <= 0) {
                iflag[0] = -3;
                throw new ExceptionWithIflag(iflag[0], "Improper input parameters  (n or m are not positive.)");
            }

            nfun = 1;
            point = 0;
            finish = false;

            if (diagco) {
                for (int i = 0; i < n; i++) {
                    if (diag[i] <= 0) {
                        iflag[0] = -2;
                        throw new ExceptionWithIflag(iflag[0], "The " + i + "-th diagonal element of the inverse hessian approximation is not positive.");
                    }
                }
            } else {
                for (int i = 0; i < n; i++) {
                    diag[i] = 1.0;
                }
            }
            ispt = n + 2 * m;
            iypt = ispt + n * m;

            for (int i = 0; i < n; i++) {
                w[ispt + i] = -g[i] * diag[i];
            }

            gnorm = Math.sqrt(ddotOpt(n, g));
            stp1 = 1 / gnorm;
            ftol = 0.0001;

            if (iprint[0] >= 0) {
                lb1(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);
            }

            execute_entire_while_loop = true;
        }

        while (true) {
            if (execute_entire_while_loop) {
                iter++;
                info[0] = 0;
                bound = iter - 1;
                if (iter != 1) {
                    if (iter > m) {
                        bound = m;
                    }
                    ys = ddotOpt2(n, w, iypt + npt, ispt + npt);
                    if (!diagco) {
                        yy = ddotOpt2(n, w, iypt + npt, iypt + npt);
                        final double invYY = 1.0/yy;

                        for (int i = 0; i < n; i++) {
                            diag[i] = ys*invYY;
                        }
                    } else {
                        iflag[0] = 2;
                        return;
                    }
                }
            }

            if (execute_entire_while_loop || iflag[0] == 2) {
                if (iter != 1) {
                    if (diagco) {
                        for (int i = 0; i < n; i++) {
                            if (diag[i] <= 0) {
                                iflag[0] = -2;
                                throw new ExceptionWithIflag(iflag[0], "The " + i + "-th diagonal element of the inverse hessian approximation is not positive.");
                            }
                        }
                    }
                    cp = point;
                    if (point == 0) {
                        cp = m;
                    }
                    w[ n + cp - 1] = 1 / ys;

                    for (int i = 0; i < n; i++) {
                        w[i] = -g[i];
                    }

                    cp = point;

                    for (int i = 0; i < bound; i++) {
                        cp -= 1;
                        if (cp == - 1) {
                            cp = m - 1;
                        }
                        sq = ddotOpt2(n, w, ispt + cp * n, 0);
                        inmc = n + m + cp;
                        iycn = iypt + cp * n;
                        w[inmc] = w[ n + cp] * sq;
                        daxpyOpt(n, -w[inmc], w, iycn);
                    }

                    for (int i = 0; i < n; i++) {
                        w[i] *= diag[i];
                    }

                    for (int i = 0; i < bound; i++) {
                        yr = ddotOpt2(n, w, iypt + cp * n, 0);
                        beta = w[ n + cp] * yr;
                        inmc = n + m + cp;
                        beta = w[inmc] - beta;
                        iscn = ispt + cp * n;
                        daxpyOpt(n, beta, w, iscn);
                        cp += 1;
                        if (cp == m) {
                            cp = 0;
                        }
                    }
                    System.arraycopy(w, 0, w, ispt + point * n, n);
                }

                nfev[0] = 0;
                stp[0] = 1;
                if (iter == 1) {
                    stp[0] = stp1;
                }
                
                System.arraycopy(g, 0, w, 0, n);
            }

            mcsrch.mcsrch(n, x, f, g, w, ispt + point * n, stp, ftol, xtol, maxfev, info, nfev, diag);

            if (info[0] == -1) {
                iflag[0] = 1;
                return;
            } else if(info[0] == 2){
                iflag[0] = -1;
                if(!beReallySilent) System.err.println("WARNING: Line search failed, interval of uncertainty is too small. This is an internal issue.");
                info[0] = -1;
                return;
            } else if(info[0] == 3){
                iflag[0] = -1;
                if(!beReallySilent) System.err.println("WARNING: Line search exceeded maximum number of function evaluations. Increase if necessary.");
                return;
            } else if(info[0] == 4){
                if(!beReallySilent) System.err.println("WARNING: LBFGS: Line search step is at the lower bound. Continuing with new gradient and energy.");
                iflag[0] = 1; // taming stuff
                info[0] = -1;
                return;
            } else if(info[0] == 5){
                if(!beReallySilent) System.err.println("WARNING: LBFGS: Line search step is at the upper bound. Continuing with new gradient and energy.");
                iflag[0] = 1; // taming stuff
                info[0] = -1;
                return;
            } else if(info[0] == 6){
                iflag[0] = -1;
                if(!beReallySilent) System.err.println("WARNING: Line search rounding errors prevent further progress. Tolerances too small?");
                info[0] = -1;
                return;
            }
            
            nfun += nfev[0];
            npt = point * n;

            for (int i = 0; i < n; i++) {
                w[ ispt + npt + i] *= stp[0];
                w[ iypt + npt + i] = g[i] - w[i];
            }

            point++;
            if (point == m) {
                point = 0;
            }

            gnorm = Math.sqrt(ddotOpt(n, g));
            xnorm = Math.sqrt(ddotOpt(n, x));
            xnorm = Math.max(1.0, xnorm);

            if (gnorm / xnorm <= eps) {
                finish = true;
            }

            if (iprint[ 1 - 1] >= 0) {
                lb1(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);
            }

            // Cache the current solution vector. Due to the spaghetti-like
            // nature of this code, it's not possible to quit here and return;
            // we need to go back to the top of the loop, and eventually call
            // mcsrch one more time -- but that will modify the solution vector.
            // So we need to keep a copy of the solution vector as it was at
            // the completion (info[0]==1) of the most recent line search.

            System.arraycopy(x, 0, solution_cache, 0, n);

            if (finish) {
                iflag[0] = 0;
                return;
            }

            execute_entire_while_loop = true;		// from now on, execute whole loop
        }
    }

    /** Print debugging and status messages for <code>lbfgs</code>.
     * Depending on the parameter <code>iprint</code>, this can include 
     * number of function evaluations, current function value, etc.
     * The messages are output to <code>System.err</code>.
     *
     * @param iprint Specifies output generated by <code>lbfgs</code>.<p>
     *		<code>iprint[0]</code> specifies the frequency of the output:
     *		<ul>
     *		<li> <code>iprint[0] &lt; 0</code>: no output is generated,
     *		<li> <code>iprint[0] = 0</code>: output only at first and last iteration,
     *		<li> <code>iprint[0] &gt; 0</code>: output every <code>iprint[0]</code> iterations.
     *		</ul><p>
     *
     *		<code>iprint[1]</code> specifies the type of output generated:
     *		<ul>
     *		<li> <code>iprint[1] = 0</code>: iteration count, number of function 
     *			evaluations, function value, norm of the gradient, and steplength,
     *		<li> <code>iprint[1] = 1</code>: same as <code>iprint[1]=0</code>, plus vector of
     *			variables and  gradient vector at the initial point,
     *		<li> <code>iprint[1] = 2</code>: same as <code>iprint[1]=1</code>, plus vector of
     *			variables,
     *		<li> <code>iprint[1] = 3</code>: same as <code>iprint[1]=2</code>, plus gradient vector.
     *		</ul>
     * @param iter Number of iterations so far.
     * @param nfun Number of function evaluations so far.
     * @param gnorm Norm of gradient at current solution <code>x</code>.
     * @param n Number of free parameters.
     * @param m Number of corrections kept.
     * @param x Current solution.
     * @param f Function value at current solution.
     * @param g Gradient at current solution <code>x</code>.
     * @param stp Current stepsize.
     * @param finish Whether this method should print the ``we're done'' message.
     */
    private static void lb1(int[] iprint, int iter, int nfun, double gnorm, int n, int m, double[] x, double f, double[] g, double[] stp, boolean finish) {

        if (iter == 0) {
            System.err.println("*************************************************");
            System.err.println("  n = " + n + "   number of corrections = " + m + "\n       initial values");
            System.err.println(" f =  " + f + "   gnorm =  " + gnorm);
            if (iprint[ 2 - 1] >= 1) {
                System.err.print(" vector x =");
                for (int i = 0; i < n; i++) {
                    System.err.print("  " + x[i]);
                }
                System.err.println("");

                System.err.print(" gradient vector g =");
                for (int i = 0; i < n; i++) {
                    System.err.print("  " + g[i]);
                }
                System.err.println("");
            }
            System.err.println("*************************************************");
            System.err.println("\ti\tnfn\tfunc\tgnorm\tsteplength");
        } else {
            if ((iprint[0] == 0) && (iter != 1 && !finish)) {
                return;
            }
            if (iprint[0] != 0) {
                if ((iter - 1) % iprint[0] == 0 || finish) {
                    if (iprint[1] > 1 && iter > 1) {
                        System.err.println("\ti\tnfn\tfunc\tgnorm\tsteplength");
                    }
                    System.err.println("\t" + iter + "\t" + nfun + "\t" + f + "\t" + gnorm + "\t" + stp[0]);
                } else {
                    return;
                }
            } else {
                if (iprint[1] > 1 && finish) {
                    System.err.println("\ti\tnfn\tfunc\tgnorm\tsteplength");
                }
                System.err.println("\t" + iter + "\t" + nfun + "\t" + f + "\t" + gnorm + "\t" + stp[0]);
            }
            if (iprint[1] == 2 || iprint[1] == 3) {
                if (finish) {
                    System.err.print(" final point x =");
                } else {
                    System.err.print(" vector x =  ");
                }
                for (int i = 0; i < n; i++) {
                    System.err.print("  " + x[i]);
                }
                System.err.println("");
                if (iprint[1] == 3) {
                    System.err.print(" gradient vector g =");
                    for (int i = 0; i < n; i++) {
                        System.err.print("  " + g[i]);
                    }
                    System.err.println("");
                }
            }
            if (finish) {
                System.err.println(" The minimization terminated without detecting errors. iflag = 0");
            }
        }
    }

    /** Compute the sum of a vector times a scalar plus another vector.
     * Adapted from the subroutine <code>daxpy</code> in <code>lbfgs.f</code>.
     * There could well be faster ways to carry out this operation; this
     * code is a straight translation from the Fortran.
     * changed this a bit, not much though (JMD)
     */
    private static void daxpyOpt(final int n, final double da, final double[] dx,
            final int ix0) {

	//assumptions: iy0 == 0 && incx == 1 && incy == 1 && dx == dy
        
        if(n <= 0 || da == 0) return;
        
        final int m = n % 4;
        if (m != 0) {
            for (int i = 0; i < m; i++) dx[i] += da * dx[ix0 + i];

            if (n < 4) return;
        }

        int t2 = ix0+m;
        for (int i = m; i < n; i++) {
            dx[i] += da * dx[t2];
            t2++;
        }
    }

    /** Compute the dot product of two vectors.
     * Adapted from the subroutine <code>ddot</code> in <code>lbfgs.f</code>.
     * There could well be faster ways to carry out this operation; this
     * code is a straight translation from the Fortran.
     * changed this a bit, not much though (JMD)
     */
    private static double ddotOpt(final int n, final double[] dx) {

        // assumptions: ix0 ==  0 && incx == 1 && iy0 == 0 && incy == 1 && dx == dy

        if (n <= 0) return 0.0;

        double d = 0.0;
        for(int i = 0; i < n; i++){
            d += dx[i]*dx[i];
        }
        
        return d;
    }
    
    /** Compute the dot product of two vectors.
     * Adapted from the subroutine <code>ddot</code> in <code>lbfgs.f</code>.
     * There could well be faster ways to carry out this operation; this
     * code is a straight translation from the Fortran.
     * changed this a bit, not much though (JMD)
     */
    private static double ddotOpt2(final int n, final double[] dx, final int ix0,
            final int iy0) {

        // assumptions: dx == dy && incx == 1 && incy == 1
        
        if (n <= 0) return 0.0;

        double d = 0.0;
        for(int i = 0; i < n; i++){
            d += dx[ix0 + i] * dx[iy0 + i];
        }
        
        return d;
        
/*        double dtemp = 0.0;
        final int m = n % 5;
        if (m != 0) {
            for (int i = 0; i < m; i++) {
                dtemp = dtemp + dx[ix0 + i] * dx[iy0 + i];
            }
            if (n < 5) {
                return dtemp;
            }
        }

        for (int i = m; i < n; i += 5) {
            dtemp = dtemp + dx[ix0 + i] * dx[iy0 + i] + dx[ix0 + i + 1] * dx[iy0 + i + 1] + dx[ix0 + i + 2] * dx[iy0 + i + 2] + dx[ix0 + i + 3] * dx[iy0 + i + 3] + dx[ix0 + i + 4] * dx[iy0 + i + 4];
        }

        return dtemp;*/
    }

    private static class Mcsrch implements Serializable {
    
        private static final long serialVersionUID = (long) 20111117;

        private final int infoc[] = new int[1];
        private final double[] dgx = new double[1], dgxm = new double[1], dgy = new double[1],
                fx = new double[1], fxm = new double[1], fy = new double[1], fym = new double[1],
                stx = new double[1], sty = new double[1], dgym = new double[1];
        private double dg = 0, dgm = 0, dginit = 0, dgtest = 0, finit = 0, ftest1 = 0,
                fm = 0, p5 = 0, p66 = 0, stmin = 0, stmax = 0, width = 0, width1 = 0, xtrapf = 0;
        private final boolean[] brackt = new boolean[1];
        private boolean stage1 = false;
        
        private final double gtol;
        private final double stepmin;
        private final double stepmax;
        
        Mcsrch(final double gtoler, final double stepmini, final double stepmaxi){
            this.gtol = gtoler;
            this.stepmin = stepmini;
            this.stepmax = stepmaxi;
        }
        
        Mcsrch(final Mcsrch orig){
            this.gtol = orig.gtol;
            this.stepmin = orig.stepmin;
            this.stepmax = orig.stepmax;
        }

        double sqr(final double x) {
            return x * x;
        }

        double max3(final double x, final double y, final double z) {
            return x < y ? (y < z ? z : y) : (x < z ? z : x);
        }

        /** Minimize a function along a search direction. This code is
         * a Java translation of the function <code>MCSRCH</code> from
         * <code>lbfgs.f</code>, which in turn is a slight modification of
         * the subroutine <code>CSRCH</code> of More' and Thuente.
         * The changes are to allow reverse communication, and do not affect
         * the performance of the routine. This function, in turn, calls
         * <code>mcstep</code>.<p>
         *
         * The Java translation was effected mostly mechanically, with some
         * manual clean-up; in particular, array indices start at 0 instead of 1.
         * Most of the comments from the Fortran code have been pasted in here
         * as well.<p>
         *
         * The purpose of <code>mcsrch</code> is to find a step which satisfies
         * a sufficient decrease condition and a curvature condition.<p>
         *
         * At each stage this function updates an interval of uncertainty with
         * endpoints <code>stx</code> and <code>sty</code>. The interval of
         * uncertainty is initially chosen so that it contains a
         * minimizer of the modified function
         * <pre>
         *      f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
         * </pre>
         * If a step is obtained for which the modified function
         * has a nonpositive function value and nonnegative derivative,
         * then the interval of uncertainty is chosen so that it
         * contains a minimizer of <code>f(x+stp*s)</code>.<p>
         *
         * The algorithm is designed to find a step which satisfies
         * the sufficient decrease condition
         * <pre>
         *       f(x+stp*s) &lt;= f(X) + ftol*stp*(gradf(x)'s),
         * </pre>
         * and the curvature condition
         * <pre>
         *       abs(gradf(x+stp*s)'s)) &lt;= gtol*abs(gradf(x)'s).
         * </pre>
         * If <code>ftol</code> is less than <code>gtol</code> and if, for example,
         * the function is bounded below, then there is always a step which
         * satisfies both conditions. If no step can be found which satisfies both
         * conditions, then the algorithm usually stops when rounding
         * errors prevent further progress. In this case <code>stp</code> only
         * satisfies the sufficient decrease condition.<p>
         *
         * @author Original Fortran version by Jorge J. More' and David J. Thuente
         *	  as part of the Minpack project, June 1983, Argonne National 
         *   Laboratory. Java translation by Robert Dodier, August 1997.
         *
         * @param n The number of variables.
         *
         * @param x On entry this contains the base point for the line search.
         *		On exit it contains <code>x + stp*s</code>.
         *
         * @param f On entry this contains the value of the objective function
         *		at <code>x</code>. On exit it contains the value of the objective
         *		function at <code>x + stp*s</code>.
         *
         * @param g On entry this contains the gradient of the objective function
         *		at <code>x</code>. On exit it contains the gradient at
         *		<code>x + stp*s</code>.
         *
         *	@param s The search direction.
         *
         * @param stp On entry this contains an initial estimate of a satifactory
         *		step length. On exit <code>stp</code> contains the final estimate.
         *
         *	@param ftol Tolerance for the sufficient decrease condition.
         *
         * @param xtol Termination occurs when the relative width of the interval
         *		of uncertainty is at most <code>xtol</code>.
         *
         *	@param maxfev Termination occurs when the number of evaluations of
         *		the objective function is at least <code>maxfev</code> by the end
         *		of an iteration.
         *
         *	@param info This is an output variable, which can have these values:
         *		<ul>
         *		<li><code>info = 0</code> Improper input parameters.
         *		<li><code>info = -1</code> A return is made to compute the function and gradient.
         *		<li><code>info = 1</code> The sufficient decrease condition and
         *			the directional derivative condition hold.
         *		<li><code>info = 2</code> Relative width of the interval of uncertainty
         *			is at most <code>xtol</code>.
         *		<li><code>info = 3</code> Number of function evaluations has reached <code>maxfev</code>.
         *		<li><code>info = 4</code> The step is at the lower bound <code>stpmin</code>.
         *		<li><code>info = 5</code> The step is at the upper bound <code>stpmax</code>.
         *		<li><code>info = 6</code> Rounding errors prevent further progress.
         *			There may not be a step which satisfies the
         *			sufficient decrease and curvature conditions.
         *			Tolerances may be too small.
         *		</ul>
         *
         *	@param nfev On exit, this is set to the number of function evaluations.
         *
         *	@param wa Temporary storage array, of length <code>n</code>.
         */
        public void mcsrch(int n, double[] x, double f, double[] g, double[] s, int is0, double[] stp, double ftol, double xtol, int maxfev, int[] info, int[] nfev, double[] wa) {
            p5 = 0.5;
            p66 = 0.66;
            xtrapf = 4;

            if (info[0] != - 1) {
                infoc[0] = 1;
                if (n <= 0 || stp[0] <= 0 || ftol < 0 || gtol < 0 || xtol < 0 || stepmin < 0 || stepmax < stepmin || maxfev <= 0) {
                    return;
                }

                // Compute the initial gradient in the search direction
                // and check that s is a descent direction.

                dginit = 0;

                for (int j = 0; j < n; j++) {
                    dginit += g[j] * s[is0 + j];
                }


                if (dginit >= 0) {
                    System.out.println("The search direction is not a descent direction.");
                    return;
                }

                brackt[0] = false;
                stage1 = true;
                nfev[0] = 0;
                finit = f;
                dgtest = ftol * dginit;
                width = stepmax - stepmin;
                width1 = width / p5;
                System.arraycopy(x, 0, wa, 0, n);

                // The variables stx, fx, dgx contain the values of the step,
                // function, and directional derivative at the best step.
                // The variables sty, fy, dgy contain the value of the step,
                // function, and derivative at the other endpoint of
                // the interval of uncertainty.
                // The variables stp, f, dg contain the values of the step,
                // function, and derivative at the current step.

                stx[0] = 0;
                fx[0] = finit;
                dgx[0] = dginit;
                sty[0] = 0;
                fy[0] = finit;
                dgy[0] = dginit;
            }

            while (true) {
                if (info[0] != -1) {
                    // Set the minimum and maximum steps to correspond
                    // to the present interval of uncertainty.

                    if (brackt[0]) {
                        stmin = Math.min(stx[0], sty[0]);
                        stmax = Math.max(stx[0], sty[0]);
                    } else {
                        stmin = stx[0];
                        stmax = stp[0] + xtrapf * (stp[0] - stx[0]);
                    }

                    // Force the step to be within the bounds stpmax and stpmin.

                    stp[0] = Math.max(stp[0], stepmin);
                    stp[0] = Math.min(stp[0], stepmax);

                    // If an unusual termination is to occur then let
                    // stp be the lowest point obtained so far.

                    if ((brackt[0] && (stp[0] <= stmin || stp[0] >= stmax)) || nfev[0] >= maxfev - 1 || infoc[0] == 0 || (brackt[0] && stmax - stmin <= xtol * stmax)) {
                        stp[0] = stx[0];
                    }

                    // Evaluate the function and gradient at stp
                    // and compute the directional derivative.
                    // We return to main program to obtain F and G.

                    for (int j = 0; j < n; j++) {
                        x[j] = wa[j] + stp[0] * s[is0 + j];
                    }

                    info[0] = -1;
                    return;
                }

                info[0] = 0;
                nfev[0] = nfev[0] + 1;
                dg = 0;

                for (int j = 0; j < n; j++) {
                    dg += g[j] * s[is0 + j];
                }

                ftest1 = finit + stp[0] * dgtest;

                // Test for convergence.

                if ((brackt[0] && (stp[0] <= stmin || stp[0] >= stmax)) || infoc[0] == 0) {
                    info[0] = 6;
                }

                if (stp[0] == stepmax && f <= ftest1 && dg <= dgtest) {
                    info[0] = 5;
                }

                if (stp[0] == stepmin && (f > ftest1 || dg >= dgtest)) {
                    info[0] = 4;
                }

                if (nfev[0] >= maxfev) {
                    info[0] = 3;
                }

                if (brackt[0] && stmax - stmin <= xtol * stmax) {
                    info[0] = 2;
                }

                if (f <= ftest1 && Math.abs(dg) <= gtol * (-dginit)) {
                    info[0] = 1;
                }

                // Check for termination.

                if (info[0] != 0) {
                    return;
                }

                // In the first stage we seek a step for which the modified
                // function has a nonpositive value and nonnegative derivative.

                if (stage1 && f <= ftest1 && dg >= Math.min(ftol, gtol) * dginit) {
                    stage1 = false;
                }

                // A modified function is used to predict the step only if
                // we have not obtained a step for which the modified
                // function has a nonpositive function value and nonnegative
                // derivative, and if a lower function value has been
                // obtained but the decrease is not sufficient.

                if (stage1 && f <= fx[0] && f > ftest1) {
                    // Define the modified function and derivative values.

                    fm = f - stp[0] * dgtest;
                    fxm[0] = fx[0] - stx[0] * dgtest;
                    fym[0] = fy[0] - sty[0] * dgtest;
                    dgm = dg - dgtest;
                    dgxm[0] = dgx[0] - dgtest;
                    dgym[0] = dgy[0] - dgtest;

                    // Call cstep to update the interval of uncertainty
                    // and to compute the new step.

                    mcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc);

                    // Reset the function and gradient values for f.

                    fx[0] = fxm[0] + stx[0] * dgtest;
                    fy[0] = fym[0] + sty[0] * dgtest;
                    dgx[0] = dgxm[0] + dgtest;
                    dgy[0] = dgym[0] + dgtest;
                } else {
                    // Call mcstep to update the interval of uncertainty
                    // and to compute the new step.

                    mcstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc);
                }

                // Force a sufficient decrease in the size of the
                // interval of uncertainty.

                if (brackt[0]) {
                    if (Math.abs(sty[0] - stx[0]) >= p66 * width1) {
                        stp[0] = stx[0] + p5 * (sty[0] - stx[0]);
                    }
                    width1 = width;
                    width = Math.abs(sty[0] - stx[0]);
                }
            }
        }

        /** The purpose of this function is to compute a safeguarded step for
         * a linesearch and to update an interval of uncertainty for
         * a minimizer of the function.<p>
         * 
         * The parameter <code>stx</code> contains the step with the least function
         * value. The parameter <code>stp</code> contains the current step. It is
         * assumed that the derivative at <code>stx</code> is negative in the
         * direction of the step. If <code>brackt[0]</code> is <code>true</code> 
         * when <code>mcstep</code> returns then a
         * minimizer has been bracketed in an interval of uncertainty
         * with endpoints <code>stx</code> and <code>sty</code>.<p>
         * 
         * Variables that must be modified by <code>mcstep</code> are 
         * implemented as 1-element arrays.
         *
         * @param stx Step at the best step obtained so far. 
         *   This variable is modified by <code>mcstep</code>.
         * @param fx Function value at the best step obtained so far. 
         *   This variable is modified by <code>mcstep</code>.
         * @param dx Derivative at the best step obtained so far. The derivative
         *   must be negative in the direction of the step, that is, <code>dx</code>
         *   and <code>stp-stx</code> must have opposite signs. 
         *   This variable is modified by <code>mcstep</code>.
         * 
         * @param sty Step at the other endpoint of the interval of uncertainty.
         *   This variable is modified by <code>mcstep</code>.
         * @param fy Function value at the other endpoint of the interval of uncertainty.
         *   This variable is modified by <code>mcstep</code>.
         * @param dy Derivative at the other endpoint of the interval of
         *   uncertainty. This variable is modified by <code>mcstep</code>.
         * 
         * @param stp Step at the current step. If <code>brackt</code> is set
         *   then on input <code>stp</code> must be between <code>stx</code>
         *   and <code>sty</code>. On output <code>stp</code> is set to the
         *   new step.
         * @param fp Function value at the current step.
         * @param dp Derivative at the current step.
         * 
         * @param brackt Tells whether a minimizer has been bracketed.
         *   If the minimizer has not been bracketed, then on input this
         *   variable must be set <code>false</code>. If the minimizer has
         *   been bracketed, then on output this variable is <code>true</code>.
         * 
         * @param stpmin Lower bound for the step.
         * @param stpmax Upper bound for the step.
         * 
         * @param info On return from <code>mcstep</code>, this is set as follows:
         *   If <code>info</code> is 1, 2, 3, or 4, then the step has been
         *   computed successfully. Otherwise <code>info</code> = 0, and this
         *   indicates improper input parameters.
         *
         * @author Jorge J. More, David J. Thuente: original Fortran version,
         *   as part of Minpack project. Argonne Nat'l Laboratory, June 1983.
         *   Robert Dodier: Java translation, August 1997.
         */
        public void mcstep(double[] stx, double[] fx, double[] dx, double[] sty, double[] fy, double[] dy, double[] stp, double fp, double dp, boolean[] brackt, double stpmin, double stpmax, int[] info) {
            boolean bound;
            double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;

            info[0] = 0;

            if ((brackt[0] && (stp[0] <= Math.min(stx[0], sty[0]) || stp[0] >= Math.max(stx[0], sty[0]))) || dx[0] * (stp[0] - stx[0]) >= 0.0 || stpmax < stpmin) {
                return;
            }

            // Determine if the derivatives have opposite sign.

            sgnd = dp * (dx[0] / Math.abs(dx[0]));

            if (fp > fx[0]) {
                // First case. A higher function value.
                // The minimum is bracketed. If the cubic step is closer
                // to stx than the quadratic step, the cubic step is taken,
                // else the average of the cubic and quadratic steps is taken.

                info[0] = 1;
                bound = true;
                theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
                s = max3(Math.abs(theta), Math.abs(dx[0]), Math.abs(dp));
                gamma = s * Math.sqrt(sqr(theta / s) - (dx[0] / s) * (dp / s));
                if (stp[0] < stx[0]) {
                    gamma = -gamma;
                }
                p = (gamma - dx[0]) + theta;
                q = ((gamma - dx[0]) + gamma) + dp;
                r = p / q;
                stpc = stx[0] + r * (stp[0] - stx[0]);
                stpq = stx[0] + ((dx[0] / ((fx[0] - fp) / (stp[0] - stx[0]) + dx[0])) / 2) * (stp[0] - stx[0]);
                if (Math.abs(stpc - stx[0]) < Math.abs(stpq - stx[0])) {
                    stpf = stpc;
                } else {
                    stpf = stpc + (stpq - stpc) / 2;
                }
                brackt[0] = true;
            } else if (sgnd < 0.0) {
                // Second case. A lower function value and derivatives of
                // opposite sign. The minimum is bracketed. If the cubic
                // step is closer to stx than the quadratic (secant) step,
                // the cubic step is taken, else the quadratic step is taken.

                info[0] = 2;
                bound = false;
                theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
                s = max3(Math.abs(theta), Math.abs(dx[0]), Math.abs(dp));
                gamma = s * Math.sqrt(sqr(theta / s) - (dx[0] / s) * (dp / s));
                if (stp[0] > stx[0]) {
                    gamma = -gamma;
                }
                p = (gamma - dp) + theta;
                q = ((gamma - dp) + gamma) + dx[0];
                r = p / q;
                stpc = stp[0] + r * (stx[0] - stp[0]);
                stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
                if (Math.abs(stpc - stp[0]) > Math.abs(stpq - stp[0])) {
                    stpf = stpc;
                } else {
                    stpf = stpq;
                }
                brackt[0] = true;
            } else if (Math.abs(dp) < Math.abs(dx[0])) {
                // Third case. A lower function value, derivatives of the
                // same sign, and the magnitude of the derivative decreases.
                // The cubic step is only used if the cubic tends to infinity
                // in the direction of the step or if the minimum of the cubic
                // is beyond stp. Otherwise the cubic step is defined to be
                // either stpmin or stpmax. The quadratic (secant) step is also
                // computed and if the minimum is bracketed then the the step
                // closest to stx is taken, else the step farthest away is taken.

                info[0] = 3;
                bound = true;
                theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
                s = max3(Math.abs(theta), Math.abs(dx[0]), Math.abs(dp));
                gamma = s * Math.sqrt(Math.max(0, sqr(theta / s) - (dx[0] / s) * (dp / s)));
                if (stp[0] > stx[0]) {
                    gamma = -gamma;
                }
                p = (gamma - dp) + theta;
                q = (gamma + (dx[0] - dp)) + gamma;
                r = p / q;
                if (r < 0.0 && gamma != 0.0) {
                    stpc = stp[0] + r * (stx[0] - stp[0]);
                } else if (stp[0] > stx[0]) {
                    stpc = stpmax;
                } else {
                    stpc = stpmin;
                }
                stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
                if (brackt[0]) {
                    if (Math.abs(stp[0] - stpc) < Math.abs(stp[0] - stpq)) {
                        stpf = stpc;
                    } else {
                        stpf = stpq;
                    }
                } else {
                    if (Math.abs(stp[0] - stpc) > Math.abs(stp[0] - stpq)) {
                        stpf = stpc;
                    } else {
                        stpf = stpq;
                    }
                }
            } else {
                // Fourth case. A lower function value, derivatives of the
                // same sign, and the magnitude of the derivative does
                // not decrease. If the minimum is not bracketed, the step
                // is either stpmin or stpmax, else the cubic step is taken.

                info[0] = 4;
                bound = false;
                if (brackt[0]) {
                    theta = 3 * (fp - fy[0]) / (sty[0] - stp[0]) + dy[0] + dp;
                    s = max3(Math.abs(theta), Math.abs(dy[0]), Math.abs(dp));
                    gamma = s * Math.sqrt(sqr(theta / s) - (dy[0] / s) * (dp / s));
                    if (stp[0] > sty[0]) {
                        gamma = -gamma;
                    }
                    p = (gamma - dp) + theta;
                    q = ((gamma - dp) + gamma) + dy[0];
                    r = p / q;
                    stpc = stp[0] + r * (sty[0] - stp[0]);
                    stpf = stpc;
                } else if (stp[0] > stx[0]) {
                    stpf = stpmax;
                } else {
                    stpf = stpmin;
                }
            }

            // Update the interval of uncertainty. This update does not
            // depend on the new step or the case analysis above.

            if (fp > fx[0]) {
                sty[0] = stp[0];
                fy[0] = fp;
                dy[0] = dp;
            } else {
                if (sgnd < 0.0) {
                    sty[0] = stx[0];
                    fy[0] = fx[0];
                    dy[0] = dx[0];
                }
                stx[0] = stp[0];
                fx[0] = fp;
                dx[0] = dp;
            }

            // Compute the new step and safeguard it.

            stpf = Math.min(stpmax, stpf);
            stpf = Math.max(stpmin, stpf);
            stp[0] = stpf;

            if (brackt[0] && bound) {
                if (sty[0] > stx[0]) {
                    stp[0] = Math.min(stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
                } else {
                    stp[0] = Math.max(stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
                }
            }
        }
    }
}
