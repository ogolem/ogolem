/**
 Copyright (c) 2019, M. Dittner, B. Hartke
               2020, D. Behrens, M. Dittner, B. Hartke
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package org.ogolem.locopt;

import java.io.Serializable;
import java.util.Arrays;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.dense.row.mult.VectorVectorMult_DDRM;

import static org.ogolem.core.Constants.ANGTOBOHR;

import org.ogolem.core.FixedValues;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericAbstractLocOpt;
import org.ogolem.generic.GenericBackend;

/**
 * Fast Inertial Relaxation Engine as a specific MD like relaxation scheme. Used
 * and benchmarked e.g. for the NEB optimization in Ref. Sheppard et al., J.
 * Chem. Phys. 2008, 158, 134106, based on the general description in: Bitzek et
 * al., Phys. Rev. Lett. 2006, 97, 170201.
 * <p>
 * We again were inspired by the (TS-)ASE version(s).
 * <p>
 * Note (NEB): Even very many different lbfgs settings did not lead to a type of
 * "blackbox" optimization scheme, i.e. severe kinks can appear. I.e. we need a
 * more direct smooth optimizer, thus FIRE. (Remark: In upper first paper, the
 * recommended optimizer is a specific lbfgs. This can indeed be tuned to
 * specific problems (finite difference first inverse diagonals), history of
 * diags and scaling here and then etc.). But overall this is not that stable
 * and can lead to erratic jumps and thus kinks in NEB optimization. For specific
 * problems, shrinking the number ot iterations and threshold etc. one can get a
 * shortcut, but not generally as "black box" optimizer.
 *
 * @author Mark Dittner
 * @author Dominik Behrens
 * @version 2020-04-27
 */
public class FIRELocOpt<E, T extends ContinuousProblem<E>> extends GenericAbstractLocOpt<E, T> {

    private final static boolean DEBUG = false;

    private static final long serialVersionUID = (long) 20200427;

    // FIRE parameters
    final private double dt;
    final private double dMaxMove;
    final private double dTMax;
    final private int nMin;
    final private double finc;
    final private double fdec;
    final private double aStart;
    final private double fa;
    final private double a;

    // conv. threshold and max iterations
    private final int iMaxIterations;
    private final double dFMax;
    private final int iNoOfMaxTrialsUntilResets;
    private final int iNoOfMaxResets;

    private final boolean doResetsToStable;
    private final boolean doResetsToBestPointSoFar;

    private final boolean historyMode;

    /**
     * Constructs a FIRE LocOpt from a pre-initialized configuration.
     *
     * @param config the FIRE configuration.
     * @param back   the Backend that is going to be used.
     */
    public FIRELocOpt(final FireConfig config, final GenericBackend<E, T> back) {
        super(back);
        this.dt = config.dt;
        this.dMaxMove = config.dMaxMove;
        this.dTMax = config.dtMax;
        this.nMin = config.nMin;
        this.finc = config.finc;
        this.fdec = config.fdec;
        this.aStart = config.astart;
        this.fa = config.fa;
        this.a = config.a;
        this.iMaxIterations = config.iMaxIterations;
        this.dFMax = config.dFMax;
        this.historyMode = back.supportsHistory();
        this.iNoOfMaxResets = config.iNoOfMaxResets;
        this.iNoOfMaxTrialsUntilResets = config.iNoOfMaxTrialsUntilResets;
        this.doResetsToStable = config.doResetToStable;
        this.doResetsToBestPointSoFar = config.doResetToBestPointSoFar;
    }

    public FIRELocOpt(final FIRELocOpt<E, T> orig) {
        super(orig);
        this.dt = orig.dt;
        this.dMaxMove = orig.dMaxMove;
        this.dTMax = orig.dTMax;
        this.nMin = orig.nMin;
        this.finc = orig.finc;
        this.fdec = orig.fdec;
        this.aStart = orig.aStart;
        this.fa = orig.fa;
        this.a = orig.a;
        this.iMaxIterations = orig.iMaxIterations;
        this.dFMax = orig.dFMax;
        this.historyMode = orig.historyMode;
        this.iNoOfMaxResets = orig.iNoOfMaxResets;
        this.iNoOfMaxTrialsUntilResets = orig.iNoOfMaxTrialsUntilResets;
        this.doResetsToStable = orig.doResetsToStable;
        this.doResetsToBestPointSoFar = orig.doResetsToBestPointSoFar;
    }

    @Override
    public FIRELocOpt<E, T> clone() {
        return new FIRELocOpt<>(this);
    }

    @Override
    public String getMyID() {
        return "Fast Internal Relaxation Engine (FIRE)";
    }

    @Override
    protected T optimize(T individual) {
        @SuppressWarnings("unchecked") final T work = (T) individual.clone();
        final int dims = back.numberOfActiveCoordinates(work);
        final double[] p = back.getActiveCoordinates(work).clone();

        final double[] gradient = new double[dims];
        double lastE = back.gradient(p, gradient, 0);

        if (DEBUG) {
            System.out.println("DEBUG: Energy in LBFGS is " + lastE);
        }

        // FIRE iterations
        boolean converged = false;

        int iNoOfResetsTaken = 0;

        Scr scr = new Scr(dims, this.dt, this.a);
        final BestPoint currentBest = historyMode ? new BestPoint() : null;
        // if resets take place, remember the point in order to investigate, if the backend really resets something (or
        // to save redundant, but often really expensive gradient evaluations)
        final double[] pointCache = this.doResetsToStable ? p.clone() : null;
        // flag: set to true if a new best point in history
        boolean didRememberANewBestPoint = false;
        boolean earlyEndBecauseOfOscillation = false;

        for (int iter = 0; iter < iMaxIterations; iter++) {
            final double convMeasure = calcConvMeasure(gradient);

            converged = convMeasure <= dFMax;
            if (DEBUG) {
                System.out.println("DEBUG: Convergence measure in FIRE is " + convMeasure);
            }

            boolean resetThisStep = false;
            if (historyMode) {
                // minus 1 as this is still the point etc. of last iteration
                boolean foundNewBest = rememberBestPoint(currentBest, p, gradient, convMeasure, lastE, iter - 1);
                if (foundNewBest && !didRememberANewBestPoint) {
                    didRememberANewBestPoint = true;
                }
                if (currentBest.iNoOfFailedTrails > iNoOfMaxTrialsUntilResets) {
                    currentBest.iNoOfFailedTrails = 0;
                    currentBest.lastResetBestQuality = convMeasure;

                    if (this.doResetsToBestPointSoFar) {
                        // an older point in history will be loaded (but the velocity(==gradient/curvature history)
                        // is reset, too (--> first steps more or less gradient descent ones, again)
                        if (didRememberANewBestPoint) {
                            // at least once between last reset and this reset a new best point was found,
                            // resetting is meaningful
                            back.historyResetToBestPoint(p, gradient);
                            didRememberANewBestPoint = false;
                        } else {
                            // no new best point was found. Without early end: same point is reloaded
                            // that was reloaded already, and velocities are also reset --> endless cycle
                            earlyEndBecauseOfOscillation = true;
                        }
                    }

                    if (this.doResetsToStable) {
                        back.resetToStable(p);
                        if (pointDidChange(pointCache, p)) {
                            lastE = back.gradient(p, gradient, iter);
                        }
                    }
                    resetThisStep = true;
                    iNoOfResetsTaken++;
                    if (DEBUG) {
                        // >> DmB Comment:
                        // This was too verbose in non-Debug-mode (its not relevant to the user anyway)
                        System.out.println("DEBUG: Too many attempts at improving target in FIRE. Reset "
                                + iNoOfResetsTaken + " taking place!");
                    }
                }
            }
            if (converged) {
                break;
            }

            if (iNoOfResetsTaken > this.iNoOfMaxResets) {
                // also break here, usually not converged until now
                break;
            }

            if (earlyEndBecauseOfOscillation) {
                System.out.println("FIRE: FIRE tackled an oscillation problem: Reloading old point with zero" +
                        "start velocities will end at the exact same state again (and again)." +
                        " Returning best point so far!");
                break;
            }

            try {
                doStep(iter, p, gradient, scr, resetThisStep);
            } catch (Exception e) {
                System.err.println("FIRE: Fatal Error in FIRE for ID " + work.getID() + ". Problem is: "
                        + e.toString());
                if (DEBUG) {
                    System.out.println("DEBUG: Convergence measure at FIRE exit: " + convMeasure + " Energy: "
                            + lastE);
                }
                break;
            }

            if (this.doResetsToStable) {
                cacheCurrentPoint(p, pointCache);
            }

            // evaluate a new gradient for this step
            lastE = back.gradient(p, gradient, iter);

            if (DEBUG) {
                System.out.println("DEBUG: Current FIRE iteration: " + iter + " Energy: " + lastE);
            }

            if (DEBUG) {
                for (int i = 0; i < dims; i++) {
                    if (p[i] != p[i]) {
                        System.err.println("DEBUG: NaN in FIRE at " + i + " iteration " + iter);
                        lastE = FixedValues.NONCONVERGEDENERGY;
                        break;
                    }
                }
            }

            if (DEBUG) {
                for (int i = 0; i < gradient.length; i++) {
                    if (gradient[i] != gradient[i]) {
                        System.err.println("DEBUG: NaN in FIRE at " + i + " iteration " + iter);

                        lastE = FixedValues.NONCONVERGEDENERGY;
                        break;
                    }
                }
            }
        }

        if (!converged) {
            final double convMeasure = calcConvMeasure(gradient);

            if (historyMode && iMaxIterations > 0 && convMeasure > currentBest.currentQualityMeasure) {
                // Comment: last check is of course redundant, usually, except for debugging
                // and "single point" evaluations using the utilities
                System.out.println("FIRE: No convergence achieved for individual " + work.getID()
                        + ". Best measure: " + currentBest.currentQualityMeasure + " of iteration: " +
                        currentBest.iteration);

                @SuppressWarnings("unchecked") final T res = (T) individual.clone();
                res.setFitness(currentBest.dEnergy);
                back.historyGetBestPoint(res);

                return res;
            } else {
                System.out.println("FIRE: No convergence achieved for individual " + work.getID());
            }
        }

        @SuppressWarnings("unchecked") final T res = (T) individual.clone();
        res.setFitness(lastE);
        back.updateActiveCoordinates(res, p);

        return res;

    }

    private double calcConvMeasure(final double[] gradient) {
        return calcConvergenceMeasureRMSG(gradient);
    }

    private void cacheCurrentPoint(final double[] p, final double[] pointCache) {
        System.arraycopy(p, 0, pointCache, 0, p.length);
    }

    private boolean pointDidChange(final double[] pointCache, final double[] p) {
        return !Arrays.equals(pointCache, p);
    }

    /**
     * Remember the best point found so far w.r.t. the quality measure (not
     * necessarily the total energy)
     *
     * @param lastE return true, if a new _best_ point was found (and remembered)
     */
    private boolean rememberBestPoint(BestPoint best, final double[] point, final double[] gradient,
                                      final double convMeasure, final double lastE, final int iteration) {
        if (convMeasure < best.currentQualityMeasure) {
            // new "global" best point found
            best.iteration = iteration;
            best.currentQualityMeasure = convMeasure;
            best.dEnergy = lastE;
            back.historyRememberPoint(point, gradient, lastE);
            best.iNoOfFailedTrails = 0; // reset
            if (DEBUG) {
                System.out.println("DEBUG: New best point in FIRE: " + lastE);
            }
            return true; // something remembered here --> true
        } else if (convMeasure < best.lastResetBestQuality) {
            // after reset, some improvements are successful, but not yet reaching the global best point above
            best.lastResetBestQuality = convMeasure;
            best.iNoOfFailedTrails = 0;
            return false;
        } else {
            // just count up a "fail" iteration
            best.iNoOfFailedTrails++; // count up
            return false;
        }
    }

    private double calcConvergenceMeasureRMSG(final double[] grad) {
        DMatrixRMaj gradMat = DMatrixRMaj.wrap(1, grad.length, grad);
        final double innerProd = VectorVectorMult_DDRM.innerProd(gradMat, gradMat) / grad.length;

        return Math.sqrt(innerProd);
    }

    /**
     * Do a FIRE step. (Actually quite easy, using Euler method instead of
     * another MD integrator). Picture: "Imagine a blind skier who is just able
     * to retard and steer looking for the fastest way to the bottom of a valley
     * in an unknown mountain range (PES)," mixing in some acceleration into the
     * direction that is steeper than the current motion (history based
     * reasonable "average" of last directions), and abruptly stopping and
     * restarting, if the new current steepest direction (force) is contrary to
     * such average motion. 1.) calc. P = F * v 2.) set v = (1-alpha)v + alpha
     * F|v| 3.1.) if P > 0 and accumulated steps since last neg. P is > nMin,
     * increase time step and decrease alpha 3.2.) if P <= 0, decrease time
     * step, freeze system and reset alpha 5.) return to MD
     *
     * @param iter     iteration number
     * @param gradient corresponding scratch object
     */
    private void doStep(final int iter, final double[] p, final double[] gradient, final Scr scr,
                        final boolean bResetThisStep) {
        DMatrixRMaj gradMat = DMatrixRMaj.wrap(1, gradient.length, gradient);
        DMatrixRMaj pMat = DMatrixRMaj.wrap(1, p.length, p);
        DMatrixRMaj veloc = scr.veloc;
        final int iNSteps = scr.iCurrStepsTaken;
        if (iter > 0) {
            final double vDot = -VectorVectorMult_DDRM.innerProd(gradMat, veloc); // -grad * v (the Power P)

            if (vDot > 0.0 && !bResetThisStep) {
                final double gradNorm = NormOps_DDRM.normF(gradMat);
                final double velocNorm = NormOps_DDRM.normF(veloc);
                DMatrixRMaj gradNormalized = scr.scrMat;
                CommonOps_DDRM.divide(gradMat, gradNorm, gradNormalized);

                // v = (1-a)*v + a * f_norm *||v||
                CommonOps_DDRM.add((1.0 - scr.a), veloc, -scr.a * velocNorm, gradNormalized, veloc);
                if (iNSteps > this.nMin) {
                    scr.dt = Math.min(scr.dt * this.finc, this.dTMax);
                    scr.a *= this.fa;
                }
                scr.iCurrStepsTaken++;
            } else {
                CommonOps_DDRM.scale(0, veloc); // all velocities set to zero
                scr.a = this.aStart;
                scr.dt *= this.fdec;
                scr.iCurrStepsTaken = 0;
                if (DEBUG) {
                    System.out.println("DEBUG: Reset velocities in FIRE at iteration " + iter);
                }
            }
        }
        DMatrixRMaj dr = scr.dr;
        CommonOps_DDRM.add(veloc, -scr.dt, gradMat, veloc); // v += dt * (-grad)

        CommonOps_DDRM.scale(scr.dt, veloc, dr); // dr = dt * v
        final double drNorm = NormOps_DDRM.fastNormF(dr); // ||dr||
        if (drNorm > this.dMaxMove) {
            // both lines together: dr = maxmove * dr/||dr||
            if (DEBUG) {
                System.out.println("DEBUG: Step in FIRE too large, rescaling!");
            }
            NormOps_DDRM.normalizeF(dr);
            CommonOps_DDRM.scale(this.dMaxMove, dr);
        }
        CommonOps_DDRM.scale(ANGTOBOHR, dr);

        CommonOps_DDRM.add(pMat, dr, pMat); // new point: p = p + dr
    }

    public FireConfig getMyConf() {
        FireConfig conf = new FireConfig();
        conf.iMaxIterations = this.iMaxIterations;
        conf.dFMax = this.dFMax;
        conf.dMaxMove = this.dMaxMove;
        // FIRE has many parameters, but most of them (benchmarked in lit.) can
        // take default values that are quite reasonable
        conf.dt = this.dt;
        conf.dtMax = this.dTMax;
        conf.nMin = this.nMin;
        conf.finc = this.finc;
        conf.fdec = this.fdec;
        conf.astart = this.aStart;
        conf.fa = this.fa;
        conf.a = this.a;

        // These are not FIRE parameters, but general history backend paramters!
        conf.iNoOfMaxResets = this.iNoOfMaxResets;
        conf.iNoOfMaxTrialsUntilResets = this.iNoOfMaxTrialsUntilResets;
        conf.doResetToStable = this.doResetsToStable;
        conf.doResetToBestPointSoFar = this.doResetsToBestPointSoFar;
        return conf;
    }

    private static class Scr implements Serializable {
        private static final long serialVersionUID = (long) 20200427;
        private final DMatrixRMaj veloc;
        private final DMatrixRMaj dr;
        private final DMatrixRMaj scrMat;
        // counter for steps taken from last freeze
        private int iCurrStepsTaken = 0;
        // variable, not just starting, integrator step length (delta time)
        private double dt;
        // variable, not just starting, alpha/a (mix)
        private double a;

        // The scratch space is made up of vectors of iProblemSize.
        public Scr(final int iProblemSize, final double startdt, final double starta) {
            this.veloc = new DMatrixRMaj(1, iProblemSize);
            this.dr = new DMatrixRMaj(1, iProblemSize);
            this.scrMat = new DMatrixRMaj(1, iProblemSize);
            this.dt = startdt;
            this.a = starta;
        }

    }

    private static class BestPoint implements Serializable {
        private static final long serialVersionUID = (long) 20200427;
        // "global" best value found so far
        double currentQualityMeasure = Double.MAX_VALUE;
        double dEnergy = Double.MAX_VALUE;
        int iteration = -1;
        int iNoOfFailedTrails = 0;
        // neg. value as this is activated after first reset, sep. best value _after_ a reset
        double lastResetBestQuality = -Double.MAX_VALUE;
    }

    /**
     * A configuration object for the FIRELocOpt. All of the parameters are initialized to their default state. Check
     * the individual fields' descriptions for more information.
     */
    public static final class FireConfig implements Serializable {
        private static final long serialVersionUID = (long) 20200427;
        /**
         * max. local optimization steps
         */
        public int iMaxIterations = 100;
        /**
         * convergence crit. for a successful opt.
         */
        public double dFMax = 1E-6;
        /**
         * max move in one step for rescaling, recommended default: 0.2
         */
        public double dMaxMove = 0.2;
        // FIRE has many parameters, but most of them (benchmarked in lit.) can
        // take default values that are quite reasonable
        /**
         * (starting) (time) step size, recommended default: 0.1
         */
        public double dt = 0.1;
        /**
         * max step (time) step size for the growing step, recommended default: 1.0
         */
        public double dtMax = 1.0;
        /**
         * number of steps without step increase, recommended default: 5
         */
        public int nMin = 5;
        /**
         * mult. factor for dt increase, recommended default: 1.1
         */
        public double finc = 1.1;
        /**
         * mult. factor for dt decrease, recommended default: 0.5
         */
        public double fdec = 0.5;
        /**
         * starting acceleration mixing factor, recommended default: 0.1
         */
        public double astart = 0.1;
        /**
         * mult. factor for acceleration mixing factor increase, recommended default: 0.99
         */
        public double fa = 0.99;
        /**
         * current acceleration mixing factor, default: 0.1
         */
        public double a = 0.1;

        // These are not FIRE parameters, but general history backend parameters!
        /**
         * max. number of resets, if too many trials lead in a wrong direction
         */
        public int iNoOfMaxResets = 1;
        /**
         * max. number of non-improving steps until doing a reset
         */
        public int iNoOfMaxTrialsUntilResets = 20;
        /**
         * whether a reset to any stable point (defined by back) should take place during optimization (with one
         * additional evaluation of that reset point)
         */
        public boolean doResetToStable = true;
        /**
         * whether a reset to the best point remembered during the optimization should take place during optimization
         * (best point and corr. gradient with new zero velocity, i.e., some anew gradient descent steps and new
         * velocity/acceleration accumulation)
         */
        public boolean doResetToBestPointSoFar = false;

        /**
         * Check whether the current config will at least produce a working optimization backend
         * (cannot test for actually good parameters of course).
         * <p>
         * This is meant to be called in the Input Readin, as it will produce Runtime Exceptions!
         */
        public void checkSanity(GenericBackend<?, ? extends ContinuousProblem<?>> back) throws RuntimeException {
            // Check if all values are within acceptable margins
            if (iMaxIterations < 0)
                throw new RuntimeException("negative maxiter in FIRE: " + iMaxIterations);
            else if (dFMax <= 0)
                throw new RuntimeException("fmax in FIRE must be positive, non-zero: " + dFMax);
            else if (dMaxMove <= 0)
                throw new RuntimeException("maxmove in FIRE must be positive, non-zero: " + dMaxMove);
            else if (dt <= 0)
                throw new RuntimeException("dt in FIRE must be positive, non-zero: " + dt);
            else if (dtMax <= 0)
                throw new RuntimeException("dtmax in FIRE must be positive, non-zero: " + dtMax);
            else if (nMin <= 0)
                throw new RuntimeException("nmin in FIRE must be positive, non-zero: " + nMin);
            else if (finc <= 0)
                throw new RuntimeException("finc in FIRE must be positive, non-zero: " + finc);
            else if (fdec <= 0)
                throw new RuntimeException("fdec in FIRE must be positive, non-zero: " + fdec);
            else if (astart <= 0)
                throw new RuntimeException("astart in FIRE must be positive, non-zero: " + astart);
            else if (fa <= 0)
                throw new RuntimeException("fa in FIRE must be positive, non-zero: " + fa);
            else if (a <= 0)
                throw new RuntimeException("a in FIRE must be positive, non-zero: " + a);
            else if (iNoOfMaxResets < 0)
                throw new RuntimeException("tryresets in FIRE must be positive: " + iNoOfMaxResets);
            else if (iNoOfMaxTrialsUntilResets < 0)
                throw new RuntimeException("trials in FIRE must be positive: " + iNoOfMaxTrialsUntilResets);
            else if (!back.supportsHistory() && doResetToBestPointSoFar)
                throw new RuntimeException("Resetting to best point in FIRE won't work without a history capable" +
                        " backend!");
        }
    }
}
