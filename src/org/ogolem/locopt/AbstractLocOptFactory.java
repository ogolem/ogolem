/**
Copyright (c) 2015, J. M. Dieterich and B. Hartke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * All advertising materials mentioning features or use of this software
      must display the following acknowledgement:

      This product includes software of the ogolem.org project developed by
      J. M. Dieterich and B. Hartke (Christian-Albrechts-University Kiel, Germany)
      and contributors.

    * Neither the name of the ogolem.org project, the University of Kiel
      nor the names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.ogolem.core.FixedValues;
import org.ogolem.generic.ContinuousProblem;
import org.ogolem.generic.GenericAbsolutePrescreeningLocOpt;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.generic.GenericNoLocOpt;

/**
 * An abstract setup for a local optimization factory.
 * @author Johannes Dieterich
 * @version 2015-04-06
 */
public abstract class AbstractLocOptFactory<E,T extends ContinuousProblem<E>> implements Serializable {
    
    private static final long serialVersionUID = (long) 20150406;
    
    private final HashMap<String,GenericBackend<E,T>> backendLookup;
    protected final double absConvThreshold;
    protected final int maxIter;
    
    protected AbstractLocOptFactory(final double absConvThreshold, final int maxIter,
            final HashMap<String,GenericBackend<E,T>> backendLookup){
        this.absConvThreshold = absConvThreshold;
        this.maxIter = maxIter;
        this.backendLookup = backendLookup;
    }
    
    public GenericLocOpt<E,T> buildLocalOpt(final String input) throws Exception {
        
        // first try to lookup specialized local optimizations
        final GenericLocOpt<E,T> special = buildSpecializedLocalOptimization(input);
        if(special != null){
            // found a specialized backend for this
            return special;
        }
        
        // now go over the generic ones
        if(input.startsWith("none:")){
            final String[] opts = tokenizeSecondLevel(input.substring(5).trim());
            
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in single point only!");
                }
            }
            
            checkSanityBackend(backend);
            
            return new GenericNoLocOpt<>(backend);
            
        } else if(input.startsWith("chained:")){
            
            final String[] opts = tokenizeFirstLevel(input.substring(8).trim());
            List<GenericLocOpt<E,T>> optims = new ArrayList<>();
            double cutoff = FixedValues.NONCONVERGEDENERGY;
            for(final String opt : opts){
                final GenericLocOpt<E,T> optim = buildLocalOpt(opt);
                optims.add(optim);
            }
            
            return new GenericChainedLocalOpt<>(optims,cutoff);
        } else if(input.startsWith("relprescreen:")){
            
            System.err.println("INFO: relative prescreening local optimization only works in shared-memory situations!");
            throw new RuntimeException("Not currently implemented: relprescreen");
            /*TODO;
            return new GenericRelativePrescreeningLocOpt();*/
        } else if(input.startsWith("absprescreen:")){
            
            final String[] split = tokenizeFirstLevel(input);
            final String[] opts = tokenizeSecondLevel(split[0]);
            
            final GenericLocOpt<E,T> secondary = buildLocalOpt(split[1]);
            
            double cutoff = FixedValues.NONCONVERGEDENERGY;
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("cutoff=")){
                    cutoff = doubleToken(opt,"cutoff=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in Apache CG!");
                }
            }
            
            checkSanityBackend(backend);
            if(secondary == null){throw new RuntimeException("Secondary local optimization in prescreening is null.");}
            
            return new GenericAbsolutePrescreeningLocOpt<>(secondary,backend,cutoff);
        } else if(input.startsWith("apachecg:")){
            
            final String[] opts = tokenizeSecondLevel(input.substring(9).trim());
            
            int noIter = maxIter;
            double absConvThresh = absConvThreshold;
            boolean useFletcher = false;
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("maxiter=")){
                    noIter = integerToken(opt,"maxiter=");
                } else if(opt.equalsIgnoreCase("usefletcher=")){
                    useFletcher = booleanToken(opt,"usefletcher=");
                } else if(opt.startsWith("absconvthresh=")){
                    absConvThresh = doubleToken(opt,"absconvthresh=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in Apache CG!");
                }
            }
            
            checkSanityConvergence(absConvThresh);
            checkSanityMaxIter(noIter);
            checkSanityBackend(backend);
            
            return new CGLocOpt<>(backend, useFletcher, absConvThresh, noIter);
        } else if(input.startsWith("lbfgs:")){
            
            final String[] opts = tokenizeSecondLevel(input.substring(6).trim());
            
            int noIter = maxIter;
            double absConvThresh = absConvThreshold;
            int noCorrections = 7;
            double gtol = 0.9;
            int lineIter = 50;
            int maxResets = 3;
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("maxiter=")){
                    noIter = integerToken(opt,"maxiter=");
                } else if(opt.startsWith("lineiter=")){
                    lineIter = integerToken(opt,"lineiter=");
                } else if(opt.startsWith("nocorrs=")){
                    noCorrections = integerToken(opt,"nocorrs=");
                } else if(opt.startsWith("gtol=")){
                    gtol = doubleToken(opt,"gtol=");
                } else if(opt.startsWith("absconvthresh=")){
                    absConvThresh = doubleToken(opt,"absconvthresh=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else if(opt.startsWith("tryresets=")){
                    maxResets = integerToken(opt,"tryresets=");
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in L-BFGS!");
                }
            }
            
            checkSanityConvergence(absConvThresh);
            checkSanityMaxIter(noIter);
            checkSanityBackend(backend);
            
            return new LBFGSLocOpt<>(backend, noCorrections, noIter, absConvThresh, gtol, lineIter, maxResets);
        } else if(input.startsWith("flemin:")){
            
            final String[] opts = tokenizeSecondLevel(input.substring(7).trim());
            
            int noIter = maxIter;
            double convThresh = absConvThreshold;
            double convThreshGrad = convThresh;
            
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("maxiter=")){
                    noIter = integerToken(opt,"maxiter=");
                } else if(opt.startsWith("convthresh=")){
                    convThresh = doubleToken(opt,"convthresh=");
                } else if(opt.startsWith("convthreshgrad=")){
                    convThreshGrad = doubleToken(opt,"convthreshgrad=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in Flemin!");
                }
            }
            
            checkSanityConvergence(convThresh);
            checkSanityConvergence(convThreshGrad);
            checkSanityMaxIter(noIter);
            checkSanityBackend(backend);
            
            return new FleminLocOpt<>(backend,convThresh,convThreshGrad,noIter);
        }  else if(input.startsWith("rnk1min:")){
            
            final String[] opts = tokenizeSecondLevel(input.substring(8).trim());
            
            int noIter = maxIter;
            double convThresh = absConvThreshold;
            double convThreshGrad = convThresh;
            double maxStep = 0.1;
            
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("maxiter=")){
                    noIter = integerToken(opt,"maxiter=");
                } else if(opt.startsWith("convthresh=")){
                    convThresh = doubleToken(opt,"convthresh=");
                } else if(opt.startsWith("convthreshgrad=")){
                    convThreshGrad = doubleToken(opt,"convthreshgrad=");
                } else if(opt.startsWith("maxstep=")){
                    maxStep = doubleToken(opt,"maxstep=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in RNK1Min!");
                }
            }
            
            checkSanityConvergence(convThresh);
            checkSanityConvergence(convThreshGrad);
            checkSanityMaxIter(noIter);
            checkSanityBackend(backend);
            
            return new RNK1MinLocOpt<>(backend,convThresh,convThreshGrad,noIter,maxStep);
        } else if(input.startsWith("fire:")) {

            final String[] opts = tokenizeSecondLevel(input.substring(5).trim());

            final FIRELocOpt.FireConfig config = new FIRELocOpt.FireConfig();
            config.iMaxIterations = maxIter;
            config.dFMax = absConvThreshold;

            GenericBackend<E, T> backend = null;
            for (final String opt : opts) {
                if (opt.startsWith("maxiter=")) {
                    config.iMaxIterations = integerToken(opt, "maxiter=");
                } else if (opt.startsWith("maxmove=")) {
                    config.dMaxMove = doubleToken(opt, "maxmove=");
                } else if (opt.startsWith("fmax=")) {
                    config.dFMax = doubleToken(opt, "fmax=");
                } else if (opt.startsWith("dt=")) {
                    config.dt = doubleToken(opt, "dt=");
                } else if (opt.startsWith("dtmax=")) {
                    config.dtMax = doubleToken(opt, "dtmax=");
                } else if (opt.startsWith("nmin=")) {
                    config.nMin = integerToken(opt, "nmin=");
                } else if (opt.startsWith("finc=")) {
                    config.finc = doubleToken(opt, "finc=");
                } else if (opt.startsWith("fdec=")) {
                    config.fdec = doubleToken(opt, "fdec=");
                } else if (opt.startsWith("astart=")) {
                    config.astart = doubleToken(opt, "astart=");
                } else if (opt.startsWith("fa=")) {
                    config.fa = doubleToken(opt, "fa=");
                } else if (opt.startsWith("a=")) {
                    config.a = doubleToken(opt, "a=");
                } else if (opt.startsWith("backend=")) {
                    backend = getBackend(opt.substring(8));
                } else if (opt.startsWith("tryresets=")) {
                    config.iNoOfMaxResets = integerToken(opt, "tryresets=");
                } else if (opt.startsWith("maxtrials=")) {
                    config.iNoOfMaxTrialsUntilResets = integerToken(opt, "maxtrials=");
                } else if (opt.startsWith("resettostable=")) {
                    config.doResetToStable = booleanToken(opt, "resettostable=");
                } else if (opt.startsWith("resettobestpoint=")) {
                    config.doResetToBestPointSoFar = booleanToken(opt, "resettobestpoint=");
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in L-BFGS!");
                }
            }

            config.checkSanity(backend);
            checkSanityBackend(backend);

            return new FIRELocOpt<>(config, backend);
        }
        
        // as a last resort...
        final GenericLocOpt<E,T> gradFreeLocOpt = buildGradFreeLocalOpt(input);
        if(gradFreeLocOpt != null){return gradFreeLocOpt;}
        
        throw new RuntimeException("No local optimization found for input: " + input + ". Tried both specialized and generic ones.");
    }
    
    protected GenericLocOpt<E,T> buildGradFreeLocalOpt(final String input) throws Exception {
        
        if(input.startsWith("bobyqa:")){
            
            final String[] opts = tokenizeSecondLevel(input.substring(7).trim());
            int noIter = maxIter;
            double initialTrust = 0.1;
            double convTrust = absConvThreshold;
            
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("maxiter=")){
                    noIter = integerToken(opt,"maxiter=");
                } else if(opt.startsWith("convtrust=")){
                    convTrust = doubleToken(opt,"convtrust=");
                } else if(opt.startsWith("inittrust=")){
                    initialTrust = doubleToken(opt,"inittrust=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in BOBYQA!");
                }
            }
            
            checkSanityConvergence(convTrust);
            checkSanityMaxIter(noIter);
            checkSanityBackend(backend);
            
            return new BOBYQALocOpt<>(backend,maxIter,initialTrust,convTrust);
        } else if(input.startsWith("newuoa:")){
            
            final String[] opts = tokenizeSecondLevel(input.substring(7).trim());
            int noIter = maxIter;
            double initialTrust = 0.1;
            double convTrust = absConvThreshold;
            
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("maxiter=")){
                    noIter = integerToken(opt,"maxiter=");
                } else if(opt.startsWith("convtrust=")){
                    convTrust = doubleToken(opt,"convtrust=");
                } else if(opt.startsWith("inittrust=")){
                    initialTrust = doubleToken(opt,"inittrust=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in BOBYQA!");
                }
            }
            
            checkSanityConvergence(convTrust);
            checkSanityMaxIter(noIter);
            checkSanityBackend(backend);
            
            return new NEWUOALocOpt<>(backend,initialTrust,convTrust,maxIter);
        } else if(input.startsWith("praxis:")){
            
            final String[] opts = tokenizeSecondLevel(input.substring(7).trim());
            
            int noIter = maxIter;
            double convThresh = absConvThreshold;
            double maxScaling = 0.1;
            double maxStep = 0.1;
            
            GenericBackend<E,T> backend = null;
            for(final String opt : opts){
                if(opt.startsWith("maxiter=")){
                    noIter = integerToken(opt,"maxiter=");
                } else if(opt.startsWith("convthresh=")){
                    convThresh = doubleToken(opt,"convthresh=");
                } else if(opt.startsWith("maxscaling=")){
                    maxScaling = doubleToken(opt,"maxscaling=");
                } else if(opt.startsWith("maxstep=")){
                    maxStep = doubleToken(opt,"maxstep=");
                } else if(opt.startsWith("backend=")){
                    backend = getBackend(opt.substring(8));
                } else {
                    throw new RuntimeException("Unknown option " + opt + " in RNK1Min!");
                }
            }
            
            checkSanityConvergence(convThresh);
            checkSanityMaxIter(noIter);
            checkSanityBackend(backend);
            
            return new PraxisLocOpt<>(backend,convThresh,maxStep, maxScaling, noIter);
        }
        
        return null;
    }
    
    protected abstract GenericLocOpt<E,T> buildSpecializedLocalOptimization(final String input) throws Exception;
    
    protected GenericBackend<E,T> getBackend(final String backend) throws Exception {
        
        if(backendLookup != null){
            final GenericBackend<E,T> ret = backendLookup.get(backend);
            if(ret != null){
                return ret;
            }
        }
        
        final GenericBackend<E,T> ret = parseBackend(backend);
        if(ret == null){
            throw new RuntimeException("No backend could be either found from either the tagged functions or through directly parsing " + backend);
        }
        
        return ret;
    }
    
    public abstract GenericBackend<E,T> parseBackend(final String backend) throws Exception;
    
    protected void checkSanityConvergence(final double convThreshold) throws Exception {
        if(convThreshold <= 0.0){
            throw new RuntimeException("Convergence threshold must be positive. Current is: " + convThreshold);
        }
    }
    
    protected void checkSanityMaxIter(final int maxIter) throws Exception {
        if(maxIter < 0){
            throw new RuntimeException("Maximum number of local optimization iterations must be at least zero. Current is: " + maxIter);
        }
    }
    
    protected void checkSanityBackend(final GenericBackend<E,T> backend) throws Exception {
        if(backend == null){
            throw new RuntimeException("Backend must be non-null.");
        }
    }
    
    /*
     * HELPER FUNCTIONS FOR OUR GRAMMAR
     */
    protected final String[] tokenizeFirstLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\|");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final String[] tokenizeSecondLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\;");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final String[] tokenizeThirdLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\,");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final String[] tokenizeFourthLevel(final String s){
        
        if(s.trim().isEmpty()){
            return new String[0];
        }
        
        final String[] tokens = s.trim().split("\\/");
        for(int i = 0; i < tokens.length; i++){
            tokens[i] = tokens[i].trim();
        }
        
        return tokens;
    }
    
    protected final boolean booleanToken(final String token, final String key){
        return Boolean.parseBoolean(token.substring(key.length()).trim());
    }
    
    protected final int integerToken(final String token, final String key){
        return Integer.parseInt(token.substring(key.length()).trim());
    }
    
    protected final double doubleToken(final String token, final String key){
        return Double.parseDouble(token.substring(key.length()).trim());
    }
}
