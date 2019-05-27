package contrib.newuoa;

/**
 * Interface defining what a method to be used with NEWUOA must handle.
 * @author Jonas Feldt
 * @author Johannes Dieterich
 * @version 2011-11-19
 */
public interface NEWUOAMethod {
    
    /**
     * Computes the function value at a defined point.
     * @param n Dimensionality of the problem.
     * @param point The point.
     * @return The function value at this point
     */
    double computeObjectiveValue(final int n, final double[] point);
}
