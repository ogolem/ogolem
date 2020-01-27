/**
Copyright (c) 2019, J. M. Dieterich and B. Hartke
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
package org.ogolem.microbenchmarks;

import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.Input;

/**
 * A library of cartesian coordinates which may be used for benchmarks.
 * @author Johannes Dieterich
 * @version 2019-12-29
 */
class CartesianCoordinatesLibrary {
    
    private static final String KANAMP2 = "69\n" +
"Coordinates from ORCA-job kana_mp2_tzvpp\n" +
"  C       0.539877     -0.718820     -1.345908\n" +
"  C       1.842787     -1.478957     -1.530449\n" +
"  C       3.003270     -0.488603     -1.510472\n" +
"  C       0.520974      0.155488     -0.096000\n" +
"  C       1.704985      1.110411     -0.081472\n" +
"  C       3.027505      0.347587     -0.235315\n" +
"  O       1.719178      1.813495      1.165317\n" +
"  O      -0.576126     -1.607260     -1.193129\n" +
"  C      -1.384640     -1.852469     -2.311705\n" +
"  C      -2.845542     -1.727565     -1.884469\n" +
"  C      -2.871644     -4.174946     -1.418625\n" +
"  C      -1.422689     -4.199451     -1.925754\n" +
"  O      -1.163813     -3.136386     -2.854913\n" +
"  C      -3.218924     -2.800245     -0.874414\n" +
"  C       1.821754      3.216921      1.054069\n" +
"  C       2.278243      3.739900      2.414517\n" +
"  C      -0.117131      4.233862      2.911889\n" +
"  C      -0.454931      3.595108      1.568432\n" +
"  O       0.618991      3.813289      0.644922\n" +
"  H       0.377796     -0.083677     -2.223230\n" +
"  H       1.947077     -2.167573     -0.680658\n" +
"  H       2.915176      0.166303     -2.383549\n" +
"  H       3.952245     -1.019037     -1.602156\n" +
"  H       0.591132     -0.493752      0.785206\n" +
"  H       1.583934      1.828559     -0.899560\n" +
"  H       3.111165     -0.321982      0.626414\n" +
"  H      -1.145144     -1.154838     -3.115863\n" +
"  H      -3.476617     -1.848726     -2.767286\n" +
"  H      -3.543839     -4.385811     -2.261002\n" +
"  H      -0.757996     -4.101916     -1.061167\n" +
"  H      -2.648738     -2.640471      0.045348\n" +
"  H       2.551242      4.788436      2.302612\n" +
"  H       0.029449      5.311738      2.771937\n" +
"  H      -0.612944      2.523360      1.702426\n" +
"  C       1.170146      3.632421      3.454577\n" +
"  H       0.969540      2.566741      3.632935\n" +
"  H       2.546252      3.481220      0.282559\n" +
"  O      -0.684751      0.903106     -0.059763\n" +
"  H      -1.401992      0.287446     -0.271964\n" +
"  N       1.808287     -2.159880     -2.825613\n" +
"  N       4.219460      1.192617     -0.232515\n" +
"  H       1.038236     -2.819328     -2.859768\n" +
"  H       2.662185     -2.687837     -2.960334\n" +
"  H       4.204490      1.815862     -1.032183\n" +
"  H       4.238452      1.777282      0.595218\n" +
"  C      -1.674443      4.197253      0.888122\n" +
"  O      -2.870040      3.976929      1.644354\n" +
"  H      -1.757296      3.776317     -0.114729\n" +
"  H      -1.564464      5.277304      0.805936\n" +
"  H      -3.051073      3.032644      1.622854\n" +
"  O       3.462164      3.067196      2.832428\n" +
"  O      -1.125289      3.986346      3.885086\n" +
"  H       3.235980      2.130169      2.837461\n" +
"  H      -1.972400      4.177883      3.464049\n" +
"  N       1.576099      4.356615      4.655856\n" +
"  H       2.418760      3.937451      5.029964\n" +
"  H       0.851417      4.259236      5.356683\n" +
"  O      -3.093196     -0.459714     -1.266388\n" +
"  H      -3.313690      0.178328     -1.950336\n" +
"  O      -3.049544     -5.108992     -0.374336\n" +
"  H      -2.676875     -5.942309     -0.703736\n" +
"  O      -4.612543     -2.763547     -0.618480\n" +
"  H      -4.801840     -1.892409     -0.257028\n" +
"  C      -1.090438     -5.468702     -2.702923\n" +
"  H      -0.090679     -5.357245     -3.132550\n" +
"  H      -1.795442     -5.555125     -3.530071\n" +
"  N      -1.233797     -6.657263     -1.848478\n" +
"  H      -1.218792     -7.495069     -2.415286\n" +
"  H      -0.448448     -6.726317     -1.212296";
    
    private static final String KANAPM3 = "69\n" +
" \n" +
"  C    0.609720698    -0.798206372    -1.459959438  \n" +
"  C    1.994419034    -1.485478547    -1.507660218  \n" +
"  C    3.090108203    -0.417904626    -1.534996023  \n" +
"  C    0.505135008     0.172107564    -0.240760977  \n" +
"  C    1.684712378     1.186533579    -0.205179979  \n" +
"  C    3.040251301     0.428504198    -0.263034288  \n" +
"  O    1.596711385     1.861219429     1.054955515  \n" +
"  O   -0.386758726    -1.810299294    -1.254143909  \n" +
"  C   -1.352750176    -1.927762211    -2.275998499  \n" +
"  C   -2.806494310    -1.741055772    -1.740561033  \n" +
"  C   -2.946246949    -4.238799565    -1.486112891  \n" +
"  C   -1.466521247    -4.308352913    -1.956029225  \n" +
"  O   -1.130415814    -3.208718819    -2.807834861  \n" +
"  C   -3.194982103    -2.885431948    -0.768859778  \n" +
"  C    1.784029118     3.263577749     1.036103863  \n" +
"  C    2.269276475     3.717402010     2.453066198  \n" +
"  C   -0.129395399     4.315097201     2.959250657  \n" +
"  C   -0.496147048     3.731877987     1.565464092  \n" +
"  O    0.580339095     3.865045105     0.639001703  \n" +
"  H    0.416240678    -0.243654076    -2.415276783  \n" +
"  H    2.114678182    -2.111838261    -0.583258158  \n" +
"  H    2.980358725     0.223448474    -2.434714196  \n" +
"  H    4.087254370    -0.896065474    -1.627422949  \n" +
"  H    0.493175342    -0.420277568     0.708975650  \n" +
"  H    1.595552770     1.899264875    -1.063191414  \n" +
"  H    3.123294088    -0.239384323     0.631215660  \n" +
"  H   -1.143388153    -1.257538822    -3.145248365  \n" +
"  H   -3.497788325    -1.744440614    -2.620894507  \n" +
"  H   -3.637855647    -4.335220551    -2.360452456  \n" +
"  H   -0.786171037    -4.300672525    -1.068202782  \n" +
"  H   -2.599040214    -2.822639064     0.175365486  \n" +
"  H    2.636988787     4.768410637     2.377412302  \n" +
"  H    0.045182904     5.417945524     2.881519764  \n" +
"  H   -0.761071613     2.648455058     1.660167103  \n" +
"  C    1.137369302     3.616903250     3.508243243  \n" +
"  H    0.902569459     2.535135290     3.698569606  \n" +
"  H    2.483433801     3.606046016     0.236881459  \n" +
"  O   -0.674549880     0.942341826    -0.320894162  \n" +
"  H   -1.402950703     0.336001072    -0.453785121  \n" +
"  N    2.109623490    -2.312398315    -2.742815652  \n" +
"  N    4.234877299     1.306355355    -0.229108126  \n" +
"  H    1.410016772    -3.024837037    -2.739187356  \n" +
"  H    3.012119395    -2.739488086    -2.763001761  \n" +
"  H    4.186050636     1.976209144    -0.969298432  \n" +
"  H    4.278802085     1.785064689     0.646466915  \n" +
"  C   -1.680378281     4.509902399     0.947548006  \n" +
"  O   -2.831226892     4.365229305     1.747269493  \n" +
"  H   -1.871502598     4.192344451    -0.094800971  \n" +
"  H   -1.507628045     5.601041258     0.946792440  \n" +
"  H   -3.165613006     3.487906750     1.601153074  \n" +
"  O    3.432548074     3.027823627     2.832684435  \n" +
"  O   -1.182420703     4.083764011     3.875479897  \n" +
"  H    3.202295645     2.114360501     2.959033593  \n" +
"  H   -1.983694369     4.350633127     3.420588993  \n" +
"  N    1.587694565     4.311715926     4.746858885  \n" +
"  H    2.361793775     3.815358470     5.135321376  \n" +
"  H    0.842690663     4.335304998     5.409959977  \n" +
"  O   -2.873983226    -0.462869826    -1.143092379  \n" +
"  H   -3.714646419    -0.407495309    -0.706026778  \n" +
"  O   -3.248955053    -5.266102400    -0.567632880  \n" +
"  H   -2.809879118    -6.056741232    -0.904342699  \n" +
"  O   -4.562001508    -2.693628392    -0.458782780  \n" +
"  H   -4.843572759    -3.460781380     0.021467258  \n" +
"  C   -1.216365541    -5.570682126    -2.806043080  \n" +
"  H   -0.195080392    -5.542244351    -3.242180674  \n" +
"  H   -1.916599729    -5.600056551    -3.665991196  \n" +
"  N   -1.448984306    -6.771075073    -1.961248900  \n" +
"  H   -1.602332717    -7.563988181    -2.547228956  \n" +
"  H   -0.647145512    -6.939997801    -1.389454007  ";
    
    static CartesianCoordinates getKanamycinAPM3Opt(){
        
        try {
            final String[] fileCont = KANAMP2.split("\n");
            final CartesianCoordinates kana = Input.parseCartesFromFileData(fileCont, 1, new int[]{69}, new short[69], new float[69]);
            
            return kana;
        } catch(Exception e){
            System.err.println("Failure to create kanamycin A MP2 cartesians. This should never happen!");
        }
        
        return null;
    }
    
    static CartesianCoordinates getKanamycinAMP2Opt(){
        
        try {
            final String[] fileCont = KANAPM3.split("\n");
            final CartesianCoordinates kana = Input.parseCartesFromFileData(fileCont, 1, new int[]{69}, new short[69], new float[69]);
            
            return kana;
        } catch(Exception e){
            System.err.println("Failure to create kanamycin A PM3 cartesians. This should never happen!");
        }
        
        return null;
    }
}
