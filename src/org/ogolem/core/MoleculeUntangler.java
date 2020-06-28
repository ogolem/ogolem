/**
 Copyright (c) 2020, D. Behrens, J. M. Dieterich, B. Hartke
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
package org.ogolem.core;

/**
 * An attempt at a simplistic approach of untangling molecules. We calculate per-atom distances, but then apply
 * gradients molecule-wise, so that molecules will simply translate as one entity.
 *
 * The force-field-like approach employs a simple gaussian function:
 *     V (r) = exp(1/(0.5*sigma*summedRadii)*x^2)
 * with the summed vdW radii of the respective atoms and user-controlled additional width parameter sigma.
 *
 * @author Dominik Behrens
 * @version 2020-06-22
 */
public class MoleculeUntangler implements CartesianFullBackend {

    private static final long serialVersionUID = (long) 20200622;

    private final double sigma;
    private final boolean moleculeWise;

    /**
     * Constructs an untangler backend.
     *
     * @param sigma Width of the gaussian part of the potential function.
     * @param moleculeWise Whether to apply gradients to entire colliding molecules, or just to individual atoms.
     */
    public MoleculeUntangler(double sigma, boolean moleculeWise) {
        assert(sigma > 0);
        this.sigma = sigma;
        this.moleculeWise = moleculeWise;
    }

    public MoleculeUntangler(MoleculeUntangler orig) {
        this.sigma = orig.sigma;
        this.moleculeWise = orig.moleculeWise;
    }

    @Override
    public CartesianFullBackend clone() {
        return new MoleculeUntangler(this);
    }

    @Override
    public String getMethodID() {
        return "Molecule Untangler (Sigma: " + sigma + ")";
    }

    @Override
    public double energyCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtomTypes, short[] atomNos,
                                    int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
                                    short[] spins, BondInfo bonds, final boolean hasRigidEnv) {
        return calculateAllTerms(iNoOfAtoms, xyz1D, atsPerMol, atomNos, null);
    }

    @Override
    public void gradientCalculation(long lID, int iIteration, double[] xyz1D, String[] saAtomTypes, short[] atomNos,
                                    int[] atsPerMol, double[] energyparts, int iNoOfAtoms, float[] faCharges,
                                    short[] spins, BondInfo bonds, Gradient grad, final boolean hasRigidEnv) {
        calculateAllTerms(iNoOfAtoms, xyz1D, atsPerMol, atomNos, grad);
    }

    private double calculateAllTerms(int iNoOfAtoms, double[] xyz1D, int[] atsPerMol, short[] atomNos, Gradient grad) {
        double[][] totGrad = null;
        double[][] derivCache = null;
        if (grad != null) {
            // In my experience, it's safer to just clean up the incoming gradient, in case it was cached by the
            // backend.
            grad.zeroGradient();
            totGrad = grad.getTotalGradient();
            derivCache = new double[2][3];
        }
        double totalEnergy = 0;
        int[][] molCache = initMolCache(atsPerMol);
        double[] radiiCache = initRadiiCache(atomNos);
        double[] pos1Cache = new double[3];
        double[] pos2Cache = new double[3];
        // This (somewhat bloated) loop goes through the upper triangle of atom-atom interactions, but molecule-wise;
        // That way it can tell which molecules it is dealing with right now (especially important for the gradients)
        for (int firstMol = 0; firstMol < molCache.length - 1; firstMol++) {
            double molEnergy = 0;
            for (int firstAtm = molCache[firstMol][0]; firstAtm < molCache[firstMol][1]; firstAtm++) {
                for (int secondMol = firstMol + 1; secondMol < molCache.length; secondMol++) {
                    for (int secondAtm = molCache[secondMol][0]; secondAtm < molCache[secondMol][1]; secondAtm++) {
                        double summedRadius = radiiCache[firstAtm] + radiiCache[secondAtm];
                        pos1Cache[0] = xyz1D[firstAtm];
                        pos1Cache[1] = xyz1D[iNoOfAtoms + firstAtm];
                        pos1Cache[2] = xyz1D[2 * iNoOfAtoms + firstAtm];
                        pos2Cache[0] = xyz1D[secondAtm];
                        pos2Cache[1] = xyz1D[iNoOfAtoms + secondAtm];
                        pos2Cache[2] = xyz1D[2 * iNoOfAtoms + secondAtm];
                        if (moleculeWise)
                            molEnergy += repulsivePotential(summedRadius, pos1Cache, pos2Cache) * (
                                    atsPerMol[firstMol] + atsPerMol[secondMol]);
                        else
                            molEnergy += repulsivePotential(summedRadius, pos1Cache, pos2Cache);
                        if (totGrad != null) {
                            repulsivePotentialDerivative(derivCache, summedRadius, pos1Cache, pos2Cache);
                            if (moleculeWise) {
                                // Add the derivatives to the whole respective molecules; can't really merge these
                                // for loops as we need two different index ranges :(
                                for (int l = molCache[firstMol][0]; l < molCache[firstMol][1]; l++) {
                                    totGrad[0][l] += derivCache[0][0];
                                    totGrad[1][l] += derivCache[0][1];
                                    totGrad[2][l] += derivCache[0][2];
                                }
                                for (int l = molCache[secondMol][0]; l < molCache[secondMol][1]; l++) {
                                    totGrad[0][l] += derivCache[1][0];
                                    totGrad[1][l] += derivCache[1][1];
                                    totGrad[2][l] += derivCache[1][2];
                                }
                            } else {
                                // Otherwise, we only add the gradient to the respective atoms, as a normal FF would.
                                totGrad[0][firstAtm] += derivCache[0][0];
                                totGrad[1][firstAtm] += derivCache[0][1];
                                totGrad[2][firstAtm] += derivCache[0][2];
                                totGrad[0][secondAtm] += derivCache[1][0];
                                totGrad[1][secondAtm] += derivCache[1][1];
                                totGrad[2][secondAtm] += derivCache[1][2];
                            }
                        }
                    }
                }
            }
            totalEnergy += molEnergy;
        }
        if (totGrad != null)
            grad.setTotalEnergy(totalEnergy);
        return totalEnergy;
    }

    /**
     * Potential energy term.
     */
    private double repulsivePotential(double summedRadius, double[] p1, double[] p2) {
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];
        double dz = p2[2] - p1[2];
        double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (dist < 1E-5) {
            // Just make sure we're consistent with the gradient function (see below)
            dist += 1E-5;
        }
        double fc = 1 / (0.5 * summedRadius * sigma);
        return Math.exp(-fc * dist * dist);
    }

    /**
     * Populates v (at least 2x3) with derivatives of the potential term.
     */
    private void repulsivePotentialDerivative(double[][] v, double summedRadius, double[] p1, double[] p2) {
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];
        double dz = p2[2] - p1[2];
        double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (dist < 1E-5) {
            // In the unlikely case that ALL atoms are exactly at 0,0,0, the gradients will be zero and the
            // untangler will do nothing. But instead of checking for this case, pretend the atoms are just a little
            // further apart than they actually are whenever two atoms are too close.
            dx += 1E-5;
            dist += 1E-5;
        }
        double fc = 1 / (0.5 * summedRadius * sigma);
        double pw = -2 * fc * Math.exp(-fc * dist * dist);
        v[0][0] = -dx * pw;
        v[0][1] = -dy * pw;
        v[0][2] = -dz * pw;
        v[1][0] = dx * pw;
        v[1][1] = dy * pw;
        v[1][2] = dz * pw;
    }

    /**
     * Small helper method to identify the starting and ending indices of each molecule's atoms.
     *
     * @param atsPerMol the atsPerMol array passed for the calculation (usually by the GenericBackend).
     * @return helper array with one entry per molecule containing lower and upper index of the molecule's atoms.
     */
    private int[][] initMolCache(int[] atsPerMol) {
        int[][] ret = new int[atsPerMol.length][2];
        int offset = 0;
        for (int i = 0; i < atsPerMol.length; i++) {
            ret[i][0] = offset;
            ret[i][1] = offset + atsPerMol[i];
            offset += atsPerMol[i];
        }
        return ret;
    }

    /**
     * Precache atom radii per index into a double array.
     */
    private double[] initRadiiCache(short[] atomNos) {
        double[] ret = new double[atomNos.length];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = AtomicProperties.giveRadius(atomNos[i]);
        }
        return ret;
    }
}
