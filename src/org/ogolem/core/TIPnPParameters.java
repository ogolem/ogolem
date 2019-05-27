/**
Copyright (c) 2014, J. M. Dieterich
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
package org.ogolem.core;

import java.io.Serializable;

/**
 * A collection of value objects containing TIP4P parameters. Numerical values
 * from wiki's Water_model page, in turn taken from Jorgensens paper of similar
 * title.
 * @author Johannes Dieterich
 * @version 2014-12-29
 */
public class TIPnPParameters implements Serializable {
    
    private static final long serialVersionUID = (long) 20150104;
    
    public static interface TIP3PParameters extends Serializable {
        
        /**
         * the O-H distance
         * @return the distance in a.u.
         */
        double getOH();
        
        /**
         * the H-O-H angle
         * @return angle in degree
         */
        double getHOH();
        
        /**
         * the LJ-A parameter (12er term) on the O
         * @return the parameter in a.u.
         */
        double getLJA();
        
        /**
         * the LJ-B parameter (6er term) on the O
         * @return the parameter in a.u.
         */
        double getLJB();
        
        /**
         * the charge parameter on O
         * @return the charge in partial charges of the electron
         */
        double getChargeO();
        
        /**
         * the charge parameter on H
         * @return the charge in partial charges of the electron
         */
        double getChargeH();
    }
    
    public static interface TIP4PParameters extends Serializable {
        
        /**
         * the O-H distance
         * @return the distance in a.u.
         */
        double getOH();
        
        /**
         * the H-O-H angle
         * @return angle in degree
         */
        double getHOH();
        
        /**
         * the O-M parameter
         * @return the O-M distance in bohr
         */
        double getOMDist();
        
        /**
         * the LJ-A parameter (12er term) on the O
         * @return the parameter in a.u.
         */
        double getLJA();
        
        /**
         * the LJ-B parameter (6er term) on the O
         * @return the parameter in a.u.
         */
        double getLJB();
        
        /**
         * the charge parameter on M
         * @return the charge in partial charges of the electron
         */
        double getChargeM();
        
        /**
         * the charge parameter on H
         * @return the charge in partial charges of the electron
         */
        double getChargeH();
    }
    
    public static class StandardTIP3PParameters implements TIP3PParameters {
        private static final long serialVersionUID = (long) 20141229;

        @Override
        public double getOH() {
            return 0.9572*Constants.ANGTOBOHR;
        }

        @Override
        public double getHOH() {
            return 104.52;
        }
        
        @Override
        public double getLJA() {
            return (582000)*Math.pow(Constants.ANGTOBOHR,12)*Constants.KCALTOHARTREE;
        }

        @Override
        public double getLJB() {
            return (595)*Math.pow(Constants.ANGTOBOHR,6)*Constants.KCALTOHARTREE;
        }

        @Override
        public double getChargeO() {
            return -0.834;
        }

        @Override
        public double getChargeH() {
            return 0.417;
        }        
    }
    
    public static class StandardTIP4PParameters implements TIP4PParameters {
        private static final long serialVersionUID = (long) 20141229;
        
        @Override
        public double getOH() {
            return 0.9572*Constants.ANGTOBOHR;
        }

        @Override
        public double getHOH() {
            return 104.52;
        }
        
        @Override
        public double getOMDist() {
            return 0.15*Constants.ANGTOBOHR;
        }

        @Override
        public double getLJA() {
            return (600000)*Math.pow(Constants.ANGTOBOHR,12)*Constants.KCALTOHARTREE;
        }

        @Override
        public double getLJB() {
            return (610)*Math.pow(Constants.ANGTOBOHR,6)*Constants.KCALTOHARTREE;
        }

        @Override
        public double getChargeM() {
            return -1.04;
        }

        @Override
        public double getChargeH() {
            return 0.52;
        }
    }
    
    public static class TIP4P_2005 implements TIP4PParameters {
        private static final long serialVersionUID = (long) 20150104;
        
        @Override
        public double getOH() {
            return 0.9572*Constants.ANGTOBOHR;
        }

        @Override
        public double getHOH() {
            return 104.52;
        }
        
        @Override
        public double getOMDist() {
            return 0.1546*Constants.ANGTOBOHR;
        }

        @Override
        public double getLJA() {
            return (731300)*Math.pow(Constants.ANGTOBOHR,12)*Constants.KCALTOHARTREE;
        }

        @Override
        public double getLJB() {
            return (736)*Math.pow(Constants.ANGTOBOHR,6)*Constants.KCALTOHARTREE;
        }

        @Override
        public double getChargeM() {
            return -1.1128;
        }

        @Override
        public double getChargeH() {
            return 0.5564;
        }
    }
}
