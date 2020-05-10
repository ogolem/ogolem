/**
Copyright (c) 2020,  M. Dittner, B. Hartke
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

package org.ogolem.adaptive;

import java.util.List;
import java.util.Random;

import org.ogolem.generic.GenericMutation;
import org.ogolem.random.Lottery;
import org.ogolem.random.RandomUtils;


/**
 * Gaussian-based fluctuations about the current input-value, within constraints.
 *
 * @author Mark Dittner
 * @version 2020-04-29
 */
class UnaryGaussianMutation implements GenericMutation<Double, AdaptiveParameters> {

    private static final long serialVersionUID = 20140718L;

    private final Lottery r;
    private final SubMode subModus;
    private final Mode mode;
    private final double[] upper;
    private final double[] lower;
    private final double gaussShaper;

    public UnaryGaussianMutation(final Mode modus, final SubMode subModus, final double gaussShaper, final double[] low, final double[] up) {
        this.gaussShaper = gaussShaper;
        this.lower = low;
        this.r = Lottery.getInstance();
        this.upper = up;
        this.subModus = subModus;
        this.mode = modus;
    }

    public UnaryGaussianMutation(final UnaryGaussianMutation orig) {
        this.r = orig.r;
        this.lower = orig.lower.clone();
        this.upper = orig.upper.clone();
        this.gaussShaper = orig.gaussShaper;
        this.subModus = orig.subModus;
        this.mode = orig.mode;
    }

    /**
     * Computes a new double value for the current gene based on the given subMode state of this operator
     */
    private static double sampleValue(final SubMode m, final double gaussShaper, final double min, final double max, final double value) {

        switch (m) {
            case CURRENT:
                return RandomUtils.gaussDoubleAroundVal(min, max, gaussShaper, value);
            case BOUND:
                return RandomUtils.gaussDouble(min, max, gaussShaper);
            default:
                throw new RuntimeException("This is a bug: Inconsistency of internal enum values and its functionality. " +
                        "Contact the authors, please.");
        }
    }

    @Override
    public UnaryGaussianMutation clone() {
        return new UnaryGaussianMutation(this);
    }

    @Override
    public String getMyID() {
        return "generic gaussian mutation"
                + "\n\t\t mode " + mode
                + "\n\t\t with submode " + subModus.toString()
                + "\n\t\t and shaper " + gaussShaper + "\n\t\t";
    }

    @Override
    public AdaptiveParameters mutate(AdaptiveParameters orig) {

        final AdaptiveParameters mut = orig.copy();

        final Double[] params = mut.getGenomeCopy();
        final int dims = params.length;

        switch (mode) {
            case ONE: {
                // find one gene to mutate
                final int loc = r.nextInt(dims);
                assert (upper[loc] >= lower[loc]);
                params[loc] = sampleValue(subModus, gaussShaper, lower[loc], upper[loc], params[loc]);
                break;
            }
            case MULTI: {
                // find a list of genes to mutate
                final int no = r.nextInt(dims);
                final List<Integer> spots = RandomUtils.listOfPoints(no, dims);
                for (final int loc : spots) {
                    assert (upper[loc] >= lower[loc]);
                    params[loc] = sampleValue(subModus, gaussShaper, lower[loc], upper[loc], params[loc]);
                }
                break;
            }
            case ALL: {
                // simply mutate every single gene
                for (int loc = 0; loc < dims; loc++) {
                    assert (upper[loc] >= lower[loc]);
                    params[loc] = sampleValue(subModus, gaussShaper, lower[loc], upper[loc], params[loc]);
                }
                break;
            }
            default:
                throw new RuntimeException("Mode '" + mode + "' unknown in bounded generic mutation. This is a bug: " +
                        "Inconsistency of internal enum values and its functionality. Contact the authors, please.");
        }

        mut.setGenome(params);

        return mut;

    }

    enum Mode {
        /**
         * Only mutate one single parameter
         */
        ONE,
        /**
         * Mutate a sub-list of the parameters present
         */
        MULTI,
        /**
         * Mutate all parameters, i.e., the full genotype
         */
        ALL
    }

    enum SubMode {
        /**
         * Sample around the current input value based on the the distribution (unary mode)
         */
        CURRENT,
        /**
         * Sample a value within the distribution and between the total bounds (nullary mode)
         */
        BOUND
    }

}
