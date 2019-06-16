/* Copyright (c) 2013-2014, J. M. Dieterich
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

      This product includes software developed at Princeton University (USA)
      by J. M. Dieterich.

    * Neither the name of Princeton University nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

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
package contrib.edu.princeton.eac;

import java.io.Serializable;
import java.util.Random;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Holds information on the grid and provides access to grid points. Works with
 * arbitrary dimensionalities of grids.
 * TODO for higher dimensions and/or increments (approx 10D), the (in my opinion very elegant) recursive scheme
 * will exhaust the default java stack size. if one really wants to use a grid solution with a
 * higher dimensionality (one should reconsider), one would need to implement a stateless grid
 * essentially, this class would be the interface (and its content could live in SimpleGrid)
 * whilst you also want a StatelessGrid that provides always only the next point for the iterator
 * or a random point w/o state, keeping the array content(!) hashsum to see if it is known.
 * this is also the caveat, as I cannot come up with a smart way, how to not keep more state and
 * still efficiently exhaust all possible points (again, for a >10D grid not that big a deal) but
 * only somehwere around 50-75% (ballpark!).
 * @author Johannes Dieterich
 * @version 2014-03-07
 */
public class Grid<T> implements Serializable, Iterable<double[]> {
    
    private static final long serialVersionUID = (long) 20130507;
    private static final Logger LOG = LoggerFactory.getLogger(Grid.class);
    private final List<double[]> gridPoints;
    private final List<T> auxDataInGrid;
    private final List<Integer> unvisited;
    private final int dim;
    private final Random r;
    
    /**
     * Construct a grid.
     * @param starts The offsets in all dimensions.
     * @param ends The ends in all dimensions.
     * @param spacing The spacing of grid points in all dimensions.
     */
    public Grid(final double[] starts, final double[] ends, final double[] spacing){
        this.gridPoints = new ArrayList<>();
        this.unvisited = new ArrayList<>();
        
        assert(starts != null);
        assert(ends != null);
        assert(spacing != null);
        assert(spacing.length == starts.length);
        assert(spacing.length == ends.length);
        
        dim = starts.length;
        
        // fill gridpoints in
        final double[] p = starts.clone();
        gridder(0, p, starts, ends, spacing);
        
        // add all the unvisted points for later usage as random offsets
        for(int i = 0; i < gridPoints.size(); i++){
            unvisited.add(i);
        }
        
        if(LOG.isDebugEnabled()){
            LOG.debug("Grid setup complete. Total number of gridpoints " + gridPoints.size());
        }
        
        this.auxDataInGrid = new ArrayList<>(gridPoints.size());
        gridPoints.forEach((_item) -> {
            auxDataInGrid.add(null);
        });
        
        assert(gridPoints.size() == unvisited.size());
        this.r = new Random();
    }
    
    private Grid(final Grid<T> orig){
        this.r = new Random();
        this.dim = orig.dim;
        this.gridPoints = new ArrayList<>();
        this.unvisited = new ArrayList<>();
        this.auxDataInGrid = new ArrayList<>();
        
        for(final double[] point : orig){
            gridPoints.add(point.clone());
        }
        
        for(int i = 0; i < orig.gridPoints.size(); i++){
            unvisited.add(i);
        }
        
        for(int i = 0; i < orig.auxDataInGrid.size(); i++){
            final T dat = orig.auxDataInGrid.get(i);
            auxDataInGrid.add(dat);
        }
        
        assert(gridPoints.size() == auxDataInGrid.size());
    }
    
    @Override
    public Grid<T> clone(){
        return new Grid<>(this);
    }
    
    /**
     * Get a random grid point which so far has not been visited.
     * @return A grid point or null if all points have been visited.
     */
    public double[] getRandomGridPoint(){
        LOG.debug("Returning a random grid point. Are there any? " + (!unvisited.isEmpty()));
        if(unvisited.isEmpty()) {return null;}
        
        final int unv = r.nextInt(unvisited.size());
        final Integer whi = unvisited.get(unv);
        unvisited.remove(whi);
        
        LOG.debug("Returning grid point " + Arrays.toString(gridPoints.get(whi)));
        return gridPoints.get(whi);
    }
    
    /**
     * Gets the list of all grid points.
     * @return A list of all grid points.
     */
    public List<double[]> getAllGridPoints(){
        LOG.debug("Returning all grid points. Number: " + gridPoints.size());
        return gridPoints;
    }
    
    /**
     * Gets the list of so far unvisited grid points.
     * @return A list of integer IDs (these are internal and correspond to the list of grid points)
     */
    public List<Integer> getAllUnvisted(){
        LOG.debug("Returning all unvisted. Currently: " + unvisited.size());
        return unvisited;
    }
    
    /**
     * Resets the unvisted grid points, allowing to go over all grid points
     * again in a random fashion.
     */
    public void resetUnvisited(){
        LOG.debug("Resetting unvisited grid points.");
        // clear out and add all points again
        unvisited.clear();
        for(int i = 0; i < gridPoints.size(); i++){
            unvisited.add(i);
        }
    }
    
    /**
     * Returns an iterator over all grid points. Will ignore (and clear) the list
     * of unvisited points. Do not use in combination with random sniffing.
     * @return Iterator returning one grid point at a time.
     */
    @Override
    public Iterator<double[]> iterator(){
        LOG.debug("Returning fresh iterator and clearing unvisted.");
        // empty the unvisited list since we are about to exhaust it anyways
        unvisited.clear();
        
        return gridPoints.iterator();
    }
    
    public void attachToPoint(final double[] point, final T data) throws Exception{
        
        final int index = gridPoints.indexOf(point);
        if(index < 0 || index > gridPoints.size()){
            throw new Exception("Point for ID not found in list of grid points.");
        }
        
        this.auxDataInGrid.set(index, data);
    }
    
    public T getDataForGridPoint(final double[] point) throws Exception {
        
        final int index = gridPoints.indexOf(point);
        if(index < 0 || index > gridPoints.size()){
            throw new Exception("Point for ID not found in list of grid points.");
        }
        
         return auxDataInGrid.get(index);
    }
    
    public List<T> getListOfAllAttachments(){
        return auxDataInGrid;
    }
    
    /**
     * Essentially the core part of this object. The gridder calls itself recursively in order to sniff through all
     * dimensions and provide a total enumeration of all grid points.
     * @param which which dimension
     * @param p the current point
     * @param starts the offsets
     * @param ends the ends
     * @param incrs the increments
     */
    private void gridder(final int which, final double[] p, final double[] starts, final double[] ends, final double[] incrs){
        
        if(LOG.isDebugEnabled()){
            String s = "";
            for(final double d : p){
                s += "  " + d;
            }
            LOG.debug("Working on " + which + " at point " + s);
        }
        
        if(p[which] > ends[which]){
            LOG.debug("Resetting " + which);
            System.arraycopy(starts, which, p, which, starts.length - which);
            return;
        } // we are at the end of this variable, reset
        
        // we are dealing with the last variable
        if(which == p.length-1){
            // loop through
            p[which] = starts[which];
            while(p[which] <= ends[which]){
                gridPoints.add(p.clone());
                if(LOG.isDebugEnabled()){
                    String s = "";
                    for(final double d : p){
                        s += "\t" + d;
                    }
                    LOG.debug("Adding " + s + " to grid.");
                }
                p[which] += incrs[which];
            }
            p[which] = starts[which]; // lets reset here
        } else {
            // go downwards without touching
            gridder(which+1,p,starts,ends,incrs);
            // go back with touching
            p[which] += incrs[which];
            gridder(which,p,starts,ends,incrs);
        }
    }
    
    /**
     * Dimensionality of this grid.
     * @return the dimensionality
     */
    public int gridDim(){
        return dim;
    }
    
    /**
     * Maps this into a pseudo xyz file... ;-)
     * @return An xyz-formatted representation of this grid. Null if dim != 3!
     */
    public String[] mapToPseudoXYZ(){
        
        if(!LOG.isDebugEnabled()){
            LOG.error("TRYING TO MAP GRID TO XYZ IN NON-DEMOMODE. STOPPING.");
            System.exit(443434);
        }
        
        if(dim != 3){
            LOG.error("Mapping to pseudo xyz only works for a 3D grid. This grid is of dimension " + dim);
            return null;
        }
        
        final String[] gr = new String[gridPoints.size()+2];
        gr[0] = "" + gridPoints.size();
        gr[1] = "pseudo xyz";
        for(int i = 0; i < gridPoints.size(); i++){
            final double[] pt = gridPoints.get(i);
            if(unvisited.contains(i)){
                // not yet...
                gr[i+2] = "He\t" + pt[0] + "\t" + pt[1] + "\t" + pt[2];
            } else {
                gr[i+2] = "O\t" + pt[0] + "\t" + pt[1] + "\t" + pt[2];
            }
        }
        
        return gr;
    }
}
