/**
Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke
              2010-2014, J. M. Dieterich
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

/**
 * This just keeps a couple of {@code static final}s, representing
 * basic physical constants. Some taken from http://folk.uio.no/michalj/node72.html
 * lately (a very helpful compilation!).
 * @author Johannes Dieterich
 * @version 2014-07-13
 */
public final class Constants {

    // disallow instantiation
    private Constants(){}

    /*
     * CONSTANTS
     */
    // Planck constant in Js
    public static final double H = 6.6260709544E-34;
    // Planck constant in a.u.
    public static final double H_AU = 2.0*Math.PI;
    // Dirac constant in Js
    public static final double HBAR = 1.05457162852E-34;
    // Dirac constant in a.u.
    public static final double HBAR_AU = 1.0;
    // atomic mass unit in kg
    public static final double U = 1.660538782E-27;
    // atomic mass unit in a.u.
    public static final double U_AU = 1822.8884798405547;
    // atomic radius (bohr radius) in m
    public static final double A = 5.29177210859E-11;
    // atomic radius (bohr radius) in a.u.
    public static final double A_AU = 1.0;
    // Boltzmann constant in a.u. (no derived unit, we choose k_B as a fundamental one)
    public static final double KB_AU = 1.0;
    // ideal gas constant in a.u. (derived unit: hbar^2/(K mol a0^2 me) )
    public static final double IDEALGAS_AU = 1.907101047994109E18;

    /*
     * CONVERSION FACTORS
     */
    public static final double HARTREETOKJ = 2625.49962;
    public static final double HARTREETOKCAL = 627.509391;
    public static final double KJTOHARTREE = 0.0003808799;
    public static final double KCALTOHARTREE = 0.001593602;
    public static final double BOHRTOANG = 0.529177249;
    public static final double ANGTOBOHR = 1.889725989;
    public static final double EVTOHARTREE = 0.036749325398;
    public static final double HARTREETOEV = 27.21138386;
    public static final double NMTOHARTREE = 45.563359;
    public static final double HARTREETONM = 0.021947460;
    public static final double AUTEMPTOKELVIN = 3.1577464E5;
    public static final double KELVINTOAUTEMP = 1./AUTEMPTOKELVIN;
    public static final double AUTOS = 2.418884326505E-17;
    public static final double STOAU = 1./AUTOS;
}
