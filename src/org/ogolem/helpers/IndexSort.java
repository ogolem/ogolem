package org.ogolem.helpers;

/**
 * BELOW ROUTINES ARE DIRECTLY FROM https://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-get-a-sorted-list-of-indices-of-an-array
 * WHICH IN TURN ORIGINATES FROM http://algs4.cs.princeton.edu/23quicksort/Quick.java
 * AND THE MUCH BETTER (W.R.T. AUTOBOXING) SOLUTION OVER A COMPARATOR.
 * LOGICALLY, THEY ARE EXEMPT FROM OUR BSD LICENSE AND REMAIN PUBLIC DOMAIN AS WE OWN NO COPYRIGHT!
 * @author PU and Stackoverflow
 */
public class IndexSort {
    
    private IndexSort(){};
    
    public static void quicksort(final double[] main, final int[] index) {
        quicksort(main, index, 0, index.length - 1);
    }

    // quicksort a[left] to a[right]
    private static void quicksort(final double[] a, final int[] index, final int left, final int right) {
        if (right <= left) return;
        final int i = partition(a, index, left, right);
        quicksort(a, index, left, i-1);
        quicksort(a, index, i+1, right);
    }

    // partition a[left] to a[right], assumes left < right
    private static int partition(final double[] a, final int[] index, 
        int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (less(a[++i], a[right]))      // find item on left to swap
                {}                              // a[right] acts as sentinel
            while (less(a[right], a[--j]))      // find item on right to swap
                {if (j == left) break;}         // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(a, index, i, j);               // swap two elements into place
        }
        exch(a, index, i, right);               // swap with partition element
        return i;
    }

    // is x < y ?
    private static boolean less(final double x, final double y) {
        return (x < y);
    }

    // exchange a[i] and a[j]
    private static void exch(final double[] a, final int[] index, final int i, final int j) {
        
        final double swap = a[i];
        a[i] = a[j];
        a[j] = swap;
        final int b = index[i];
        index[i] = index[j];
        index[j] = b;
    }
}
