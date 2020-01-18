/**
Copyright (c) 2017-2018, J. M. Dieterich and B. Hartke
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

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.ogolem.adaptive.genericfitness.GenericReferencePoint;
import org.ogolem.adaptive.genericfitness.PropertyCalculator;
import org.ogolem.adaptive.genericfitness.ReferenceGenericMatrixData;
import org.ogolem.adaptive.genericfitness.ReferenceGenericScalarData;
import org.ogolem.adaptive.genericfitness.ReferenceGenericTensorData;
import org.ogolem.adaptive.genericfitness.ReferenceGenericVectorData;
import org.ogolem.adaptive.genericfitness.ReferenceInputData;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.StreamGobbler;
import org.ogolem.io.InputPrimitives;
import org.ogolem.io.ManipulationPrimitives;
import org.ogolem.io.OutputPrimitives;
import org.ogolem.properties.GenericMatrixProperty;
import org.ogolem.properties.GenericScalarProperty;
import org.ogolem.properties.GenericTensorProperty;
import org.ogolem.properties.GenericVectorProperty;
import org.ogolem.properties.Property;

/**
 * An interface to wire through the generic properties.
 * @author Johannes Dieterich
 * @version 2018-01-02
 */
public class GenericCaller implements Adaptivable {

    private static final long serialVersionUID = (long) 20180102;
    
    private final int noParams;
    private final double[][] borders;
    
    public GenericCaller(final int noParams) throws Exception {
    
        assert(noParams > 0);
        
        this.noParams = noParams;
        this.borders = new double[2][noParams];
        
        final String[] dat = InputPrimitives.readFileIn("generic.borders");
        for(int i = 0; i < noParams; i++){
            final String[] line = dat[i].trim().split("\\s+");
            borders[0][i] = Double.parseDouble(line[0]);
            borders[1][i] = Double.parseDouble(line[1]);
        }
    }
    
    private GenericCaller(final GenericCaller orig){
        this.noParams = orig.noParams;
        this.borders = orig.borders;
    }
    
    @Override
    public GenericCaller clone() {
        return new GenericCaller(this);
    }

    @Override
    public double energyOfStructWithParams(CartesianCoordinates cartes, AdaptiveParameters params, int geomID, BondInfo bonds) {
        throw new UnsupportedOperationException("Not supported. Don't use generic caller for this!");
    }

    @Override
    public double gradientOfStructWithParams(CartesianCoordinates cartes, AdaptiveParameters params, int geomID, BondInfo bonds, double[] grad) {
        throw new UnsupportedOperationException("Not supported. Don't use generic caller for this!");
    }

    @Override
    public double[][] minMaxBordersForParams(AdaptiveParameters params) {
        return borders;
    }

    @Override
    public AdaptiveParameters createInitialParameterStub(ArrayList<CartesianCoordinates> refCartes, String sMethod) {
        
        final String[] atoms = new String[]{"hartkenium"};
        final int[] paramsPerAt = new int[]{noParams};
        final AdaptiveParameters paramStub = new AdaptiveParameters(noParams, -1, atoms, paramsPerAt, sMethod);
        
        // add parameter descriptions
        final String[] descr = paramStub.getAllDescriptions();
        for(int i = 0; i < descr.length; i++){
            descr[i] = "generic parameter " + i;
        }

        return paramStub;
    }

    @SuppressWarnings("unchecked")
    @Override
    public <T extends Property, V extends ReferenceInputData<T>> PropertyCalculator<T,V> getCalculatorForProperty(final T property, final V data) {
        
        // instantiate nested classes
        if(property instanceof ScalarCalculator){
            if(data instanceof ReferenceGenericScalarData){
                return (PropertyCalculator<T,V>) new ScalarCalculator();
            }
        } else if(property instanceof VectorCalculator){
            if(data instanceof ReferenceGenericVectorData){
                return (PropertyCalculator<T,V>) new VectorCalculator();
            }
        } else if(property instanceof MatrixCalculator){
            if(data instanceof ReferenceGenericMatrixData){
                return (PropertyCalculator<T,V>) new MatrixCalculator();
            }
        } else if(property instanceof TensorCalculator){
            if(data instanceof ReferenceGenericTensorData){
                return (PropertyCalculator<T,V>) new TensorCalculator();
            }
        } 
        
        System.err.println("No calculator for property " + property.name() + " with data " + data.getClass().getName() + " in GenericCaller.");
        
        return null;
    }

    @Override
    public List<? extends Property> runAllPropertyCalcs(AdaptiveParameters params, List<GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>>> referencePoints) {
        
        final List<Property> allProps = new ArrayList<>();
        for(final GenericReferencePoint<? extends Property, ? extends ReferenceInputData<?>> refPoint : referencePoints){
            final Property p = refPoint.getReferenceProperty();
            if(p instanceof GenericScalarProperty){
                final GenericScalarProperty ref = (GenericScalarProperty) p;
                final ReferenceGenericScalarData dat = (ReferenceGenericScalarData) refPoint.getReferenceInputData();
                final PropertyCalculator<GenericScalarProperty,ReferenceGenericScalarData> calc = getCalculatorForProperty(ref, dat);
            
                final GenericScalarProperty prop = calc.calculateProperty(params, dat);
                allProps.add(prop);
            } else if(p instanceof GenericVectorProperty){
                final GenericVectorProperty ref = (GenericVectorProperty) p;
                final ReferenceGenericVectorData dat = (ReferenceGenericVectorData) refPoint.getReferenceInputData();
                final PropertyCalculator<GenericVectorProperty,ReferenceGenericVectorData> calc = getCalculatorForProperty(ref, dat);
            
                final GenericVectorProperty prop = calc.calculateProperty(params, dat);
                allProps.add(prop);
            } else if(p instanceof GenericMatrixProperty){
                final GenericMatrixProperty ref = (GenericMatrixProperty) p;
                final ReferenceGenericMatrixData dat = (ReferenceGenericMatrixData) refPoint.getReferenceInputData();
                final PropertyCalculator<GenericMatrixProperty,ReferenceGenericMatrixData> calc = getCalculatorForProperty(ref, dat);
            
                final GenericMatrixProperty prop = calc.calculateProperty(params, dat);
                allProps.add(prop);
            } else if(p instanceof GenericTensorProperty){
                final GenericTensorProperty ref = (GenericTensorProperty) p;
                final ReferenceGenericTensorData dat = (ReferenceGenericTensorData) refPoint.getReferenceInputData();
                final PropertyCalculator<GenericTensorProperty,ReferenceGenericTensorData> calc = getCalculatorForProperty(ref, dat);
            
                final GenericTensorProperty prop = calc.calculateProperty(params, dat);
                allProps.add(prop);
            } else {
                // error
                throw new RuntimeException("Unknown property type " + p.printableProperty() + ".");
            }
        }
        
        return allProps;
    }
    
    class ScalarCalculator implements PropertyCalculator<GenericScalarProperty,ReferenceGenericScalarData> {

        private static final long serialVersionUID = (long) 20180102;
        
        private final String magicScript;
        
        ScalarCalculator(){
            final String tmp = System.getenv("OGO_GENERICSCALARSCRIPT");
            if(tmp == null){
                magicScript = "ogo_scalar_caller";
            } else {
                magicScript = tmp;
            }
        }
        
        private ScalarCalculator(final ScalarCalculator orig){
            this.magicScript = orig.magicScript;
        }
        
        @Override
        public ScalarCalculator clone() {
            return new ScalarCalculator(this);
        }

        @Override
        public GenericScalarProperty calculateProperty(AdaptiveParameters p, ReferenceGenericScalarData data) {
            
            final long pointID = data.getPointID();
            final int typeID = data.getTypeID();
            
            final String folder = "scalar_" + typeID + "_" + pointID + "-" + System.currentTimeMillis();
            
            // print parameter values out into a file
            final double[] params = p.getAllParamters();
            final String[] dat = new String[params.length];
            for(int i = 0; i < params.length; i++){
                dat[i] = "" + params[i];
            }
            
            try{
                OutputPrimitives.writeOut(folder + File.separator + "parameters", dat, false);
            } catch(Exception e){
                System.err.println("Failed to write parameters for generic scalar property.");
                return new GenericScalarProperty(FixedValues.NONCONVERGEDENERGY,typeID);
            }
            
            // call a script with the ID of the scalar type and the ID of the reference point
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{magicScript, "" + typeID, "" + pointID};
                final File dir = new File(folder);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    System.err.println("Failed to wait for generic scalar property.");
                    return new GenericScalarProperty(FixedValues.NONCONVERGEDENERGY,typeID);
                }

            } catch (Exception e) {
                System.err.println("Generic scalar caller has a problem. ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folder);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up generic scalar files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new GenericScalarProperty(FixedValues.NONCONVERGEDENERGY,typeID);
            }
            
            // read in the data
            double scalar = FixedValues.NONCONVERGEDENERGY;
            try{
                final String[] out = InputPrimitives.readFileIn(folder + File.separator + "scalar.out");
                scalar = Double.parseDouble(out[0].trim());
            } catch (Exception e){
                System.err.println("Failed to read generic scalar property.");
                return new GenericScalarProperty(FixedValues.NONCONVERGEDENERGY,typeID);
            }
            
            return new GenericScalarProperty(scalar,typeID);
        }

        @Override
        public GenericScalarProperty calculatePropertyGradient(AdaptiveParameters p, ReferenceGenericScalarData data, double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }        
    }
    
    class VectorCalculator implements PropertyCalculator<GenericVectorProperty,ReferenceGenericVectorData> {

        private static final long serialVersionUID = (long) 20180102;
        
        private final String magicScript;
        
        VectorCalculator(){
            final String tmp = System.getenv("OGO_GENERICVECTORSCRIPT");
            if(tmp == null){
                magicScript = "ogo_vector_caller";
            } else {
                magicScript = tmp;
            }
        }
        
        private VectorCalculator(final VectorCalculator orig){
            this.magicScript = orig.magicScript;
        }
        
        @Override
        public VectorCalculator clone() {
            return new VectorCalculator(this);
        }

        @Override
        public GenericVectorProperty calculateProperty(AdaptiveParameters p, ReferenceGenericVectorData data) {
            
            final long pointID = data.getPointID();
            final int typeID = data.getTypeID();
            
            final String folder = "vector_" + typeID + "_" + pointID + "-" + System.currentTimeMillis();
            
            // print parameter values out into a file
            final double[] params = p.getAllParamters();
            final String[] dat = new String[params.length];
            for(int i = 0; i < params.length; i++){
                dat[i] = "" + params[i];
            }
            
            try{
                OutputPrimitives.writeOut(folder + File.separator + "parameters", dat, false);
            } catch(Exception e){
                System.err.println("Failed to write parameters for generic vector property.");
                return new GenericVectorProperty(null, true, typeID);
            }
            
            // call a script with the ID of the scalar type and the ID of the reference point
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{magicScript, "" + typeID, "" + pointID};
                final File dir = new File(folder);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    System.err.println("Failed to wait for generic vector property.");
                    return new GenericVectorProperty(null, true,typeID);
                }

            } catch (Exception e) {
                System.err.println("Generic vector caller has a problem. ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folder);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up generic vector files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new GenericVectorProperty(null, true,typeID);
            }
            
            // read in the data
            double[] vector = null;
            try{
                final String[] out = InputPrimitives.readFileIn(folder + File.separator + "vector.out");
                vector = new double[out.length];
                for(int i = 0; i < out.length; i++){
                    vector[i] = Double.parseDouble(out[i].trim());
                }
            } catch (Exception e){
                System.err.println("Failed to read generic vector property.");
                return new GenericVectorProperty(null, true,typeID);
            }
            
            return new GenericVectorProperty(vector, true,typeID);
        }

        @Override
        public GenericVectorProperty calculatePropertyGradient(AdaptiveParameters p, ReferenceGenericVectorData data, double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }        
    }
    
    class MatrixCalculator implements PropertyCalculator<GenericMatrixProperty,ReferenceGenericMatrixData> {

        private static final long serialVersionUID = (long) 20180102;
        
        private final String magicScript;
        
        MatrixCalculator(){
            final String tmp = System.getenv("OGO_GENERICMATRIXSCRIPT");
            if(tmp == null){
                magicScript = "ogo_matrix_caller";
            } else {
                magicScript = tmp;
            }
        }
        
        private MatrixCalculator(final MatrixCalculator orig){
            this.magicScript = orig.magicScript;
        }
        
        @Override
        public MatrixCalculator clone() {
            return new MatrixCalculator(this);
        }

        @Override
        public GenericMatrixProperty calculateProperty(AdaptiveParameters p, ReferenceGenericMatrixData data) {
            
            final long pointID = data.getPointID();
            final int typeID = data.getTypeID();
            
            final String folder = "matrix_" + typeID + "_" + pointID + "-" + System.currentTimeMillis();
            
            // print parameter values out into a file
            final double[] params = p.getAllParamters();
            final String[] dat = new String[params.length];
            for(int i = 0; i < params.length; i++){
                dat[i] = "" + params[i];
            }
            
            try{
                OutputPrimitives.writeOut(folder + File.separator + "parameters", dat, false);
            } catch(Exception e){
                System.err.println("Failed to write parameters for generic matrix property.");
                return new GenericMatrixProperty(null, true, typeID);
            }
            
            // call a script with the ID of the scalar type and the ID of the reference point
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{magicScript, "" + typeID, "" + pointID};
                final File dir = new File(folder);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    System.err.println("Failed to wait for generic matrix property.");
                    return new GenericMatrixProperty(null, true,typeID);
                }

            } catch (Exception e) {
                System.err.println("Generic matrix caller has a problem. ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folder);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up generic matrix files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new GenericMatrixProperty(null, true,typeID);
            }
            
            // read in the data
            double[][] matrix = null;
            try{
                final String[] out = InputPrimitives.readFileIn(folder + File.separator + "matrix.out");
                matrix = new double[out.length][];
                for(int i = 0; i < out.length; i++){
                    final String[] line = out[i].trim().split("\\s+");
                    matrix[i] = new double[line.length];
                    for(int j = 0; j < out.length; j++){
                        matrix[i][j] = Double.parseDouble(line[j]);
                    }
                }
            } catch (Exception e){
                System.err.println("Failed to read generic matrix property.");
                return new GenericMatrixProperty(null, true,typeID);
            }
            
            return new GenericMatrixProperty(matrix, true,typeID);
        }

        @Override
        public GenericMatrixProperty calculatePropertyGradient(AdaptiveParameters p, ReferenceGenericMatrixData data, double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }        
    }
    
    class TensorCalculator implements PropertyCalculator<GenericTensorProperty,ReferenceGenericTensorData> {

        private static final long serialVersionUID = (long) 20180102;
        
        private final String magicScript;
        
        TensorCalculator(){
            final String tmp = System.getenv("OGO_GENERICTENSORSCRIPT");
            if(tmp == null){
                magicScript = "ogo_tensor_caller";
            } else {
                magicScript = tmp;
            }
        }
        
        private TensorCalculator(final TensorCalculator orig){
            this.magicScript = orig.magicScript;
        }
        
        @Override
        public TensorCalculator clone() {
            return new TensorCalculator(this);
        }

        @Override
        public GenericTensorProperty calculateProperty(AdaptiveParameters p, ReferenceGenericTensorData data) {
            
            final long pointID = data.getPointID();
            final int typeID = data.getTypeID();
            
            final String folder = "tensor_" + typeID + "_" + pointID + "-" + System.currentTimeMillis();
            
            // print parameter values out into a file
            final double[] params = p.getAllParamters();
            final String[] dat = new String[params.length];
            for(int i = 0; i < params.length; i++){
                dat[i] = "" + params[i];
            }
            
            try{
                OutputPrimitives.writeOut(folder + File.separator + "parameters", dat, false);
            } catch(Exception e){
                System.err.println("Failed to write parameters for generic tensor property.");
                return new GenericTensorProperty(null, true, typeID);
            }
            
            // call a script with the ID of the scalar type and the ID of the reference point
            Process proc;
            try{
                final Runtime rt = Runtime.getRuntime();
                final String[] cmd = new String[]{magicScript, "" + typeID, "" + pointID};
                final File dir = new File(folder);
                proc = rt.exec(cmd,null,dir);

                // any error message?
                final StreamGobbler errorGobbler = new
                    StreamGobbler(proc.getErrorStream(), "ERROR");

                // any output?
                final StreamGobbler outputGobbler = new
                    StreamGobbler(proc.getInputStream(), "OUTPUT");

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                if(proc.waitFor() != 0){
                    System.err.println("Failed to wait for generic tensor property.");
                    return new GenericTensorProperty(null, true,typeID);
                }

            } catch (Exception e) {
                System.err.println("Generic tensor caller has a problem. ");
                e.printStackTrace(System.err);
                try{
                    ManipulationPrimitives.remove(folder);
                } catch(Exception e1){
                    System.err.println("WARNING: Failed to clean up generic tensor files. ");
                    e1.printStackTrace(System.err);
                }
                
                return new GenericTensorProperty(null, true,typeID);
            }
            
            // read in the data
            double[][][] tensor = null;
            try{
                final String[] out = InputPrimitives.readFileIn(folder + File.separator + "tensor.out");
                
                final String[] line = out[0].trim().split("\\s+");
                
                final int xDim = Integer.parseInt(line[0]);
                final int yDim = Integer.parseInt(line[1]);
                final int zDim = Integer.parseInt(line[2]);
                
                tensor = new double[xDim][yDim][zDim];
                
                int count = 1;
                for (int x = 0; x < xDim; x++) {
                    for (int y = 0; y < yDim; y++) {
                        for (int z = 0; z < zDim; z++) {
                            tensor[x][y][z] = Double.parseDouble(out[count]);
                            count++;
                        }
                    }
                }
            } catch (Exception e){
                System.err.println("Failed to read generic tensor property.");
                return new GenericTensorProperty(null, true,typeID);
            }
            
            return new GenericTensorProperty(tensor, true,typeID);
        }

        @Override
        public GenericTensorProperty calculatePropertyGradient(AdaptiveParameters p, ReferenceGenericTensorData data, double[] grad) {
            return NumericalGradients.numericalPropertyGradient(this, p, data, grad);
        }        
    }
}
