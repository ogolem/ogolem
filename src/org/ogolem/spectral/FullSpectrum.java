/**
Copyright (c) 2014-2015, J. M. Dieterich
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
package org.ogolem.spectral;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * A full spectrum, i.e., a spectrum on an equidistant grid.
 * @author Johannes Dieterich
 * @version 2015-03-01
 */
public class FullSpectrum implements Spectrum {
    
    private static final long serialVersionUID = (long) 20141129;
    private static final boolean DEBUG = false;
    private List<Peak> peaks = null;
    private final double offset;
    private final double incr;
    private final double[] intensities;
    
    public FullSpectrum(final double start, final double incr, final double[] intensities){
        assert(start >= 0.0);
        assert(incr >= 0.0);
        this.offset = start;
        this.incr = incr;
        this.intensities = intensities;
        normalize();
    }
    
    public FullSpectrum(final FullSpectrum orig){
        this.offset = orig.offset;
        this.incr = orig.incr;
        this.intensities = orig.intensities.clone();
    }
    
    public FullSpectrum(final List<FullSpectrum> contribs, final double[] weights){
        assert(weights.length == contribs.size());
        assert(contribs.size() >= 1);
        this.offset = contribs.get(0).offset;
        this.incr = contribs.get(0).incr;
        this.intensities = new double[contribs.get(0).intensities.length];
        
        // add up
        int c = 0;
        for(final FullSpectrum spec : contribs){
            assert(spec.incr == this.incr);
            assert(spec.offset == this.offset);
            assert(spec.intensities.length == this.intensities.length);
            final double[] inte = spec.intensities;
            final double w = weights[c];
            for(int i = 0; i < inte.length; i++){
                intensities[i] += w*inte[i];
            }
            c++;
        }
        
        // normalize
        normalize();
    }
    
    public FullSpectrum(final String[] specData) throws Exception{
        
        double mySt = -1;
        double myIncr = -1;
        final List<Double> intes = new ArrayList<>();
        for(final String s : specData){
            if(s.startsWith("#")){continue;}
            final String[] sa = s.trim().split("\\s+");
            if(mySt == -1){
                mySt = Double.parseDouble(sa[0]);
                intes.add(Double.parseDouble(sa[1]));
                continue;
            }
            if(myIncr == -1){
                myIncr = Double.parseDouble(sa[0]) - mySt;
                intes.add(Double.parseDouble(sa[1]));
                continue;
            }
            intes.add(Double.parseDouble(sa[1]));
        }
        
        this.offset = mySt;
        this.incr = myIncr;
        this.intensities = new double[intes.size()];
        for(int i = 0; i < intes.size(); i++){
            intensities[i] = intes.get(i);
        }
        
        normalize();
    }
    
    @Override
    public FullSpectrum clone(){
        return new FullSpectrum(this);
    }
    
    @Override
    public final void normalize(){
        
        // find the maximal peak intensity
        double maxInt = 0.0;
        for(int i = 0; i < getIntensities().length; i++){
            maxInt = Math.max(maxInt, getIntensities()[i]);
        }
        
        if(Math.abs(maxInt) < 1e-8){
            return;
        }
        
        // normalize to that
        for(int i = 0; i < getIntensities().length; i++){
            intensities[i] /= maxInt;
        }
    }
    
    public final void normalizeToBaseline(){
        
        double minInte = Double.MAX_VALUE;
        for(final double d : intensities){
            minInte = Math.min(minInte, d);
        }
            
        for(int i = 0; i < intensities.length; i++){
            intensities[i] -= minInte;
        }
        
        normalize();
    }
    
    @Override
    public String getFormattedSpectrum(){
        
        String out = "";
        for(int i = 0; i < getIntensities().length; i++){
            out += String.format(Locale.US, "%8.1f      %7.5e\n", (getOffset()+getIncr()*i), getIntensities()[i]);
        }
        
        return out;
    }

    /**
     * @return the offset
     */
    public double getOffset() {
        return offset;
    }

    /**
     * @return the incr
     */
    public double getIncr() {
        return incr;
    }

    /**
     * @return the intensities
     */
    public double[] getIntensities() {
        return intensities;
    }

    @Override
    public FullSpectrum addUp(final List<Spectrum> spectra, final double[] weights, final double threshPeakLocEqual,
            final double coeffThresh) throws Exception {
        
        List<FullSpectrum> all = new ArrayList<>();
        for(final Spectrum spec : spectra){
            if(spec instanceof FullSpectrum){
                all.add((FullSpectrum) spec);
            } else {
                throw new Exception("All spectra must be of implementation FullSpectrum, at least this one is not.");
            }
        }
        
        return new FullSpectrum(all, weights);
    }
    
    public FullSpectrum addUp(final List<Spectrum> spectra, final double[] weights) throws Exception {
        return addUp(spectra, weights, -1.0, -1.0); // as for the full spectrum case, this anyways does not matter
    }

    @Override
    public List<Peak> getPeaks() {
        if(this.peaks != null){ return peaks;}
        
        this.normalize();
        this.peaks = findPeaks();
        return peaks;
    }
    
    @Override
    public List<Peak> getPeaksInInterval(final double start, final double end){
        
        final List<Peak> peaksInter = new ArrayList<>();
        // we assume stuff is sorted by now
        final List<Peak> allPeaks = getPeaks();
        for(final Peak p : allPeaks){
            if(p.getLocation() >= start && p.getLocation() <= end){
                peaksInter.add(p.clone());
            }
        }
        
        return peaksInter;
    }
    
    private List<Peak> findPeaks(){
        
        final List<Peak> peaksLoc = new ArrayList<>();
        
        final double thresh = SpectrumConfig.diffInPeakThreshold;
        assert(thresh >= 0.0);
        
        double loc = offset;
        double lastLoc = 0.0;
        double lastInte = 0.0;
        boolean downHill = false;
        double startInte = 0.0;
        double platStart = 0.0;
        for(int i = 0; i < intensities.length; i++){
            final double inte = intensities[i];
            if(inte >= thresh){
                if(DEBUG){System.out.println("DEBUG: Status for " + loc + " lastLoc "
                        + lastLoc + " lastInte " + lastInte + " startInte " + startInte + " platStart " + platStart
                        + " downHill " + downHill + " inte " + inte + " thresh " + thresh);}
                if((inte-lastInte) > thresh || (inte-startInte) > thresh){
                    // uphill: clear cut
                    lastLoc = loc;
                    platStart = loc;
                    downHill = false;
                    startInte = inte;
                    lastInte = inte;
                    if(DEBUG){System.out.println("DEBUG: uphill " + lastLoc + " " + lastInte);}
                } else if(Math.abs(inte-lastInte) < thresh && Math.abs(inte-startInte) < thresh){
                    // we are on a plateau
                    lastLoc = loc;
                    lastInte = inte;
                    //downHill = false; // set it to false to be sure a downhill is properly detected
                    // do NOT touch downHill. That way we can make sure that if we are on a downhill slope
                    // we do not mistakenly detect another peak by setting it to false here
                    if(DEBUG){System.out.println("DEBUG: on plateau " + lastLoc + " " + lastInte);}
                } else if((lastInte-inte) > thresh || (startInte-inte) > thresh){
                    if(!downHill){
                        // first downhill -> peak!
                        final double peakLoc = (lastLoc+platStart)/2;
                        final double peakInte = (lastInte+startInte)/2;
                        if(DEBUG){System.out.println("DEBUG: Peak located at " + peakLoc + " with intensity " + peakInte);}
                        final Peak p = new Peak(peakLoc,peakInte);
                        peaksLoc.add(p);
                        downHill = true;
                        lastLoc = loc;
                        lastInte = inte;
                        startInte = inte;
                        platStart = loc;
                    } else {
                        downHill = true;
                        lastLoc = loc;
                        lastInte = inte;
                        startInte = inte;
                        platStart = loc;
                        if(DEBUG){System.out.println("DEBUG: Going downhill " + lastInte + " and " + inte);}
                    }
                } else {
                    System.err.println("ERROR: Weird case in findPeaks()");
                    throw new Error("Weird case in findPeaks()");
                }
            }
            loc += incr;
        }
        
        return peaksLoc;
    }
    
    public FullSpectrum getSubSpectrumNormed(final double start, final double end, final boolean toBaseLine) {
        
        final List<Double> tmp = new ArrayList<>();
        double curr = offset;
        for(int i = 0; i < intensities.length; i++){
            if(curr < start){
                curr += incr;
                continue;
            } else if(curr > end){
                break;
            } else {
                tmp.add(intensities[i]);
            }
            curr += incr;
        }
        
        final double[] intes = new double[tmp.size()];
        for(int i = 0; i < tmp.size(); i++){
            intes[i] = tmp.get(i);
        }
        
        final FullSpectrum subSpec = new FullSpectrum(start,incr,intes);
        if(toBaseLine){
            subSpec.normalizeToBaseline();
        } else {
            subSpec.normalize();
        }
        
        return subSpec;
    }
}
