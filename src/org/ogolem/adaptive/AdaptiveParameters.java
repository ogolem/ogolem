/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
              2015, J. M. Dieterich and B. Hartke
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
package org.ogolem.adaptive;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import org.ogolem.generic.ContinuousProblem;

/**
 * A value object holding the parameters. Since this is the central point for a
 * good performance, some caching is build-in. It is important to note that the
 * cache is not thread-safe. Therefore in a multi-threaded environment, it is
 * absolutely required to generate a copy of this class for each thread.
 * @author Johannes Dieterich
 * @version 2015-10-28
 */
public class AdaptiveParameters extends ContinuousProblem<Double> {

    // for serialization purposes
    private static final long serialVersionUID = (long) 20150730;

    private final HashMap<String,int[]> hashPosition;
    private final String forMethod;
    private final String[] keyIDs;
    private final int[] paramsPerKey;
    private final double[] parameters;
    private final String[] paramDescriptions;

    private final HashMap<String,double[]> hashCache;

    /**
     * A constructor, just needing the the number of free parameters.
     * @param numberOfParameters the total number of parameters
     * @param id the id of this set
     * @param keys
     * @param noOfParamsPerKey
     * @param paramsForMethod
     */
    public AdaptiveParameters(final int numberOfParameters, final int id, final String[] keys,
            final int[] noOfParamsPerKey, final String paramsForMethod){
        super();
        this.parameters = new double[numberOfParameters];
        this.paramDescriptions = new String[numberOfParameters];
        this.myID = id;
        this.forMethod = paramsForMethod;
        this.hashPosition = new HashMap<>();
        this.hashCache = new HashMap<>();
        this.keyIDs = keys;
        this.paramsPerKey = noOfParamsPerKey;

        // populate the hashtables
        int iOffset = 0;
        for(int i = 0; i < keys.length; i++){
            int iAmount = noOfParamsPerKey[i];
            final int[] iaPositions = {iAmount,iOffset};
            this.hashPosition.put(keys[i], iaPositions);
            iOffset = iOffset + iAmount;
        }
        
    }

    /**
     * A copy constructor
     * @param original
     */
    public AdaptiveParameters(final AdaptiveParameters original){
        super(original);
        parameters = original.getAllParamtersCopy();
        forMethod = original.getForWhichMethod();
        hashPosition = original.getCompletePositionsCopy();
        keyIDs = original.getForWhichAtomsCopy();
        paramsPerKey = original.paramsPerKey.clone();
        // we *on purpose* empty the hashCache
        hashCache = new HashMap<>();
        paramDescriptions = original.paramDescriptions.clone();
    }

    public AdaptiveParameters(final AdaptiveParameters original, final boolean shallow) {
        super(original);
        if (shallow) {
            parameters = original.getAllParamtersCopy();
            forMethod = original.getForWhichMethod();
            keyIDs = original.getForWhichAtomsCopy();
            paramsPerKey = original.paramsPerKey.clone();
            hashPosition = original.getCompletePositionsShallowCopy();
            // we *on purpose* empty the hashCache
            hashCache = new HashMap<>();
            paramDescriptions = original.paramDescriptions;
        } else {
            parameters = original.getAllParamtersCopy();
            forMethod = original.getForWhichMethod();
            hashPosition = original.getCompletePositionsCopy();
            keyIDs = original.getForWhichAtomsCopy();
            paramsPerKey = original.paramsPerKey.clone();
            // we *on purpose* empty the hashCache
            hashCache = new HashMap<>();
            paramDescriptions = original.paramDescriptions.clone();
        }
    }

    /**
     * Create the parameter object from formatted input.
     * @param formattedParams the formatted parameters as in a .param file.
     * @param id the id of this set
     */
    public AdaptiveParameters(final String[] formattedParams, final int id){
        super();
        this.myID = id;

        String sTempForMethod = null;
        double dFitnessVal = 0.0;
        final List<Double> llParams = new LinkedList<>();
        final List<String> llParamDecrs = new LinkedList<>();
        final List<String> alForAtoms = new ArrayList<>();

        int iParamOffset = 0;

        hashCache = new HashMap<>();
        hashPosition = new HashMap<>();

        for(int i = 0; i < formattedParams.length; i++){
            String sTemp = formattedParams[i].trim();
            if(sTemp.equalsIgnoreCase("<PARAMETERSFOR>")){
                // check whether there is also the closing tag two lines later
                if(! formattedParams[i+2].trim().equalsIgnoreCase("</PARAMETERSFOR>")){
                    System.err.println("ERROR: no proper information found for what method this" +
                            " set of parameters is. Setting it to null now. Continuing.");
                    sTempForMethod = null;
                    continue;
                } else {
                    sTempForMethod = formattedParams[i+1].trim();
                    i += 2;
                    continue;
                }
            } else if(sTemp.equalsIgnoreCase("<PARAMETERFITNESS>")){
                // check whether there is also the closing tag two lines later
                if(! formattedParams[i+2].trim().equalsIgnoreCase("</PARAMETERFITNESS>")){
                    System.err.println("ERROR: no proper information found on parameter fitness" +
                            " set of parameters is. Setting it to NONCONVERGEDENERGY now. Continuing.");
                    dFitnessVal = FixedValues.NONCONVERGEDENERGY;
                    continue;
                } else {
                    dFitnessVal= Double.parseDouble(formattedParams[i+1].trim());
                    i += 2;
                    continue;
                }
            } else if(sTemp.equalsIgnoreCase("<PARAMETERS>")){
                // check for the atom tag
                String sAtom;
                int iParamAmount = 0;
                if(!formattedParams[i+1].trim().equalsIgnoreCase("<ATOM>")){
                    System.err.println("ERROR: no proper information found for which atom these " +
                            "parameters are. Continuing...");
                    continue;
                } else {
                    if(!formattedParams[i+3].trim().equalsIgnoreCase("</ATOM>")){
                        System.err.println("ERROR: no proper information found for which atom these " +
                            "parameters are. Continuing. Closing tag was: " + formattedParams[i+3].trim());
                        continue;
                    } else {
                        sAtom = formattedParams[i+2].trim();
                        i += 4;
                    }
                }

                // probe how many parameters we have
                if(!formattedParams[i].trim().equalsIgnoreCase("<VALUES>")){
                    System.err.println("ERROR: no proper information found about " +
                            "parameters. Continuing.");
                        continue;
                } else {
                    i++;

                    while(!formattedParams[i].trim().equalsIgnoreCase("</VALUES>") && i < formattedParams.length){
                        
                        final String line = formattedParams[i].trim();
                        final String[] tokens = line.split("\\s+");
                        
                        double param;
                        try{
                            param = Double.parseDouble(tokens[0]);
                        } catch(Exception e){
                            System.err.println("ERROR: Casting a parameter value failed. " + e.toString() + ". Value was " + tokens[0]);
                            throw new RuntimeException("ERROR: Casting a parameter value failed. " + e.toString() + ". Value was " + tokens[0]);
                        }
                        llParams.add(param);
                        
                        // try to parse (if existing) the description. But this shall be non-lethal if it does not work!
                        String description = null;
                        if(tokens.length > 1){
                            // we assume the rest to simply be the description
                            for(int x = 1; x < tokens.length; x++){
                                description += tokens[x];
                            }
                        }
                        llParamDecrs.add(description);

                        iParamAmount++;
                        i++;
                    }

                    // one increment to get off from the </VALUES> to </PARAMETERS> so that the next increment works again
                    i++;

                    // adding info to the hashtables
                    final int[] iaPosition = {iParamAmount,iParamOffset};
                    hashPosition.put(sAtom,iaPosition);
                    alForAtoms.add(sAtom);

                    iParamOffset += iParamAmount;

                    continue;
                }
            }

            // if we got through all tests, it must be at least worth a warning
            System.err.println("WARNING: Not all lines in the parameter file were usable. " +
                    "This might indicate a problem.");
            System.err.println("WARNING: Unusable line was " + sTemp);
        }

        forMethod = sTempForMethod;

        fitness = dFitnessVal;

        keyIDs = new String[alForAtoms.size()];
        for(int i = 0; i < alForAtoms.size(); i++){
            keyIDs[i] = alForAtoms.get(i);
        }

        paramsPerKey = new int[alForAtoms.size()];
        for(int i = 0; i < paramsPerKey.length; i++){
            final int[] info = hashPosition.get(keyIDs[i]);
            paramsPerKey[i] = info[0];
        }

        // copy the parameters
        parameters = new double[llParams.size()];
        for(int i = 0; i < llParams.size(); i++){
            parameters[i] = llParams.get(i);
        }
        
        // put the parameter descriptions in
        paramDescriptions = new String[llParamDecrs.size()];
        for(int i = 0; i < llParamDecrs.size(); i++){
            // may be null, but that is fine
            paramDescriptions[i] = llParamDecrs.get(i);
        }
    }
    
    @Override
    public AdaptiveParameters clone(){
        return new AdaptiveParameters(this);
    }

    @Override
    public void setFitness(final double fitness){
        this.fitness = fitness;
    }

    @Override
    public double getFitness(){
        return this.fitness;
    }

    public int getNumberOfParamters(){
        return parameters.length;
    }

    String getForWhichMethod(){
        return forMethod;
    }

    String getForWhichKey(final int position){

        int iCurr = 0;
        for(int i = 0; i < paramsPerKey.length; i++){
            iCurr += paramsPerKey[i];
            if(position < iCurr){
                return keyIDs[i];
            }
        }

        System.err.println("ERROR: This may be a bug in AdaptiveParameters, returning null. Contact the author(s) please.");
        return null;
    }

    String[] getForWhichAtoms(){
        return keyIDs;
    }
   
    String[] getForWhichAtomsCopy(){
        return keyIDs.clone();
    }

    public String[] getAllKeysCopy(){
        return keyIDs.clone();
    }

    private HashMap<String,int[]> getCompletePositionsCopy(){
        final HashMap<String, int[]> hashCopy = new HashMap<>();

        final Set<String> allKeys = hashPosition.keySet();
        final Iterator<String> itKeys = allKeys.iterator();

        while(itKeys.hasNext()){
            String sKey = itKeys.next();
            final int[] iaPosition = hashPosition.get(sKey).clone();
            hashCopy.put(sKey, iaPosition);
        }

        return hashCopy;
    }

    private HashMap<String,int[]> getCompletePositionsShallowCopy(){

        final HashMap<String, int[]> hashCopy = new HashMap<>(hashPosition);

        return hashCopy;
    }

    /**
     * Handle with absolute care! This returns a reference to all parameters
     * be aware that the cache stays untouched. So should you manipulate a
     * value in this array and the parameter is still cached, the change is not
     * transfered completely until you use the setAllParameters() function which
     * clears the cache.
     * @return All parameters in this object as a reference.
     */
    public double[] getAllParamters(){
        return parameters;
    }

    /**
     * Returns a copy of all parameters.
     * @return A copy of all parameters.
     */
    double[] getAllParamtersCopy(){
        return parameters.clone();
    }

    /**
     * Returns the parameters for this specific key. Please do note
     * that any change in these parameters does not get back into the total
     * genome of this individual but instead just poisons the cache.
     * @param key An unique key.
     * @return The parameters for this particular key.
     */
    public double[] getParametersForKey(final String key){

        /*
         * first try the cache ;-)
         */
        final double[] cachedParams = hashCache.get(key);
        if(cachedParams != null){
            // very good
            return cachedParams;
        }

        /*
         * was not cached, try to find out if it even exists
         */
        final int[] positions = hashPosition.get(key);

        if(positions == null){
            // catch the case that no info can be found
            return null;
        }

        final double[] params = new double[positions[0]];

        System.arraycopy(parameters,positions[1],params,0,positions[0]);

        // cache it now for the next time(s)
        hashCache.put(key, params);

        return params;
    }

    int getAmountOfParametersForKey(final String key){
        return hashPosition.get(key)[0];
    }

    @Override
    public void setGenome(final Double[] genome){
        assert(genome != null);
        assert(genome.length == parameters.length);
        for(int i = 0; i < genome.length; i++){
            parameters[i] = genome[i];
        }
        clearCache();
    }
    
    /**
     * Copies parameters in and clears the internal cache.
     * @param params the parameters to be copied in. Must be of the same length (and in the same order) as the previous parameters
     */
    public void copyParameters(final double[] params){
        assert(params != null);
        assert(parameters != null);
        assert(params.length == parameters.length);
        System.arraycopy(params, 0, parameters, 0, parameters.length);
        // empty the cache
        clearCache();
    }

    void setParameterAtPos(final double param, final int pos){
        parameters[pos] = param;
        // clear the cache since it is now useless
        clearCache();
    }
    
    String[] getAllDescriptions(){
        return paramDescriptions;
    }
    
    void setDescriptionAtPos(final String descr, final int pos){
        paramDescriptions[pos] = descr;
    }

    public String[] createFormattedOutput() {
        
        final LinkedList<String> llFormatOut = new LinkedList<>();

        // first a tag saying <PARAMETERSFOR>
        llFormatOut.add("<PARAMETERSFOR>");
        llFormatOut.add(forMethod);
        llFormatOut.add("</PARAMETERSFOR>");

        // then we put the fitness there
        llFormatOut.add("<PARAMETERFITNESS>");
        llFormatOut.add(" " + fitness);
        llFormatOut.add("</PARAMETERFITNESS>");

        // then we put all atoms into a <PARAMETER> tags, where we have <ATOM> and <VALUES>
        for(final String key : keyIDs){
            llFormatOut.add("<PARAMETERS>");
            llFormatOut.add("<ATOM>");
            llFormatOut.add(key);
            llFormatOut.add("</ATOM>");

            // the actual values
            llFormatOut.add("<VALUES>");
            final int[] iaPosition = hashPosition.get(key);
            final int iOffset = iaPosition[1];
            final int iAmount = iaPosition[0];
            for (int i = iOffset; i < iOffset + iAmount; i++) {
                String line = " " + parameters[i];
                if(paramDescriptions[i] != null){
                    // we have a description for this parameter
                    line += "\t" + paramDescriptions[i];
                }
                llFormatOut.add(line);
            }
            llFormatOut.add("</VALUES>");
            llFormatOut.add("</PARAMETERS>");
        }

        // copy the information
        final String[] saFormatOut = new String[llFormatOut.size()];
        for(int i = 0; i < llFormatOut.size(); i++){
            saFormatOut[i] = llFormatOut.get(i);
        }

        return saFormatOut;
    }

    void setRandomParams(final double[] lowerBound, final double[] upperBound){
        
        assert(lowerBound.length == upperBound.length);
        assert(lowerBound.length == parameters.length);
        
        final Random random = new Random();

        for(int i = 0; i < parameters.length; i++){
            assert(upperBound[i] >= lowerBound[i]);
            final double diff = Math.abs(upperBound[i] - lowerBound[i]);
            final double rand = random.nextDouble();
            parameters[i] = lowerBound[i] + (rand * diff);
        }
        clearCache();
    }

    List<String> getKeysStartingWith(final String start){

        final List<String> corrKeys = new ArrayList<>();

        final Set<String> allKeys = hashPosition.keySet();
        for(String key : allKeys){
            if(key.startsWith(start)){
                corrKeys.add(key);
            }
        }

        return corrKeys;
    }

    /**
     * Clears the cache.
     */
    void clearCache(){
        hashCache.clear();
    }

    /**
     * Get the start point for this key.
     * @param key The key.
     * @return The start of this key. Negative if the key is not known.
     */
    public int getStartPointForKey(final String key){
        // we can't use the cache, of course
        final int[] iaPosition = hashPosition.get(key);

        if(iaPosition == null) return -1;
        else return iaPosition[1];
    }

    /**
     * Get the start and end point for parameters with the given prefix.
     * @param prefix
     * @return The start and end point for this prefix.
     */
    int[] getStartAndEndForPrefix(final String prefix){

        final int[] startEnd = new int[2];
        int counter = 0;
        boolean startFound = false;
        for (int i = 0; i < keyIDs.length; i++) {
            if (!startFound && keyIDs[i].startsWith(prefix)) {
                startFound = true;
                startEnd[0] = counter;
            }

            if (startFound && i == keyIDs.length - 1) {
                // last entry
                startEnd[1] = parameters.length - 1;

                return startEnd;
            }

            if (startFound && !keyIDs[i].startsWith(prefix)) {
                startEnd[1] = counter - 1;

                // everything fine here
                return startEnd;
            }

            counter += paramsPerKey[i];
        }

        // if we end up here, there was somthing wrong
        return null;
    }
    
    @Override
    public double[] getGenomeAsDouble(){
        return parameters;
    }
    
    @Override
    public Double[] getGenomeCopy(){
        
        final Double[] genome = new Double[parameters.length];
        for(int i = 0; i < parameters.length; i++){
            genome[i] = parameters[i];
        }
        
        return genome;
    }
    
    /**
     * Use the actual parameter values to obtain a unique ID. A bit slower but absolutely worth it. Collisions between different objects are intentional.
     * @return a unique idea for the parameter *values*
     */
    public long getUniqueID(){
        final int hashCode = Arrays.hashCode(parameters);
        return (long) hashCode;
    }
}
