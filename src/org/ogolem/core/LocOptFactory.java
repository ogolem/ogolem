/*
Copyright (c) 2015-2020, J. M. Dieterich and B. Hartke
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

import java.io.File
import java.util.HashMap;
import org.ogolem.adaptive.AdaptiveParameters;
import org.ogolem.generic.GenericBackend;
import org.ogolem.generic.GenericLocOpt;
import org.ogolem.heat.LocalHeatPulses;
import org.ogolem.interfaces.GenericInterface;
import org.ogolem.interfaces.OrcaCaller;
import org.ogolem.interfaces.XTBCaller;
import org.ogolem.locopt.AbstractLocOptFactory;

/**
 * A local optimization factory for the geometry optimizable.
 *
 * @author Johannes Dieterich
 * @version 2020-04-25
 */
public class LocOptFactory extends AbstractLocOptFactory<Molecule, Geometry> {

  private static final long serialVersionUID = (long) 2020425;

  private final GlobalConfig config;
  private final AdaptiveParameters params;

  public LocOptFactory(
      final GlobalConfig config,
      final AdaptiveParameters params,
      final HashMap<String, GenericBackend<Molecule, Geometry>> defBackends) {
    super(config.threshLocOptGradient, config.maxIterLocOpt, defBackends);
    this.config = config;
    this.params = params;
  }

  @Override
  protected GenericLocOpt<Molecule, Geometry> buildSpecializedLocalOptimization(final String input)
      throws Exception {

    if (input.startsWith("localheat:")) {
      final String s = input.substring(10).trim();
      final int i = s.indexOf(";");
      // other data
      final String s2 = s.substring(0, i);
      // the locopt info
      final String s3 = s.substring(i + 1);
      // recursive call
      final LocOptFactory factory =
          new LocOptFactory(
              config, config.parameters, config.backendDefs); // XXX parameters may not be suitable!
      final GenericLocOpt<Molecule, Geometry> n = factory.buildLocalOpt(s3);
      if (n == null) {
        throw new RuntimeException("No Locopt found for " + s3);
      }

      final LocalHeatPulses.Configuration heatConfig = new LocalHeatPulses.Configuration();
      final String[] tokens = s2.split("\\,+");
      for (int x = 0; x < tokens.length; x++) {
        final String token = tokens[x].trim();
        if (token.startsWith("choosemode=")) {
          final String sub = token.substring(11).trim();
          switch (sub) {
            case "pick5":
              heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.PICK5;
              break;
            case "10percent":
              heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.PERCENT10;
              break;
            case "upto5":
              heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.UPTO5;
              break;
            case "upto10percent":
              heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.UPTO10PERCENT;
              break;
            case "onsphere":
              heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.ONSPHERE;
              break;
            case "insphere":
              heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.INSPHERE;
              break;
            case "incenter":
              heatConfig.chooseMode = LocalHeatPulses.Configuration.CHOOSEMODE.INCENTER;
              break;
            default:
              throw new RuntimeException("Unknown move mode " + sub + " for local heat pulses.");
          }
        } else if (token.startsWith("docddd=")) {
          heatConfig.doCDDD = Boolean.parseBoolean(token.substring(7).trim());
        } else if (token.startsWith("eqiter=")) {
          heatConfig.eqIter = Integer.parseInt(token.substring(7).trim());
        } else if (token.startsWith("iters=")) {
          heatConfig.iters = Integer.parseInt(token.substring(6).trim());
        } else if (token.startsWith("movemode=")) {
          final String sx = token.substring("movemode=".length()).trim();
          if (sx.equalsIgnoreCase("coms")) {
            heatConfig.moveMode = LocalHeatPulses.Configuration.MOVEMODE.MOVECOMS;
          } else if (sx.equalsIgnoreCase("cartesian")) {
            heatConfig.moveMode = LocalHeatPulses.Configuration.MOVEMODE.MOVECARTESIAN;
          } else {
            throw new RuntimeException("Illegal move mode " + s2 + " for local heat.");
          }
        } else if (token.startsWith("scalefac=")) {
          heatConfig.scaleFac = Double.parseDouble(token.substring(9).trim());
        } else if (token.startsWith("startampl=")) {
          heatConfig.startAmplitude = Double.parseDouble(token.substring(10).trim());
        } else if (token.startsWith("temperature=")) {
          heatConfig.temperature = Double.parseDouble(token.substring(12).trim());
        } else if (token.startsWith("usemetropolis=")) {
          heatConfig.useMetropolis = Boolean.parseBoolean(token.substring(14).trim());
        } else if (token.startsWith("starteuler=")) {
          heatConfig.startEulerStrength = Double.parseDouble(token.substring(11).trim());
        } else if (token.startsWith("sigmax=")) {
          heatConfig.sigmaX = Double.parseDouble(token.substring(7).trim());
        } else if (token.startsWith("sigmay=")) {
          heatConfig.sigmaY = Double.parseDouble(token.substring(7).trim());
        } else if (token.startsWith("sigmaz=")) {
          heatConfig.sigmaZ = Double.parseDouble(token.substring(7).trim());
        } else if (token.startsWith("sigma=")) {
          final double d = Double.parseDouble(token.substring(6).trim());
          heatConfig.sigmaX = d;
          heatConfig.sigmaY = d;
          heatConfig.sigmaZ = d;
        } else if (token.equalsIgnoreCase("verbose")) {
          heatConfig.printVerbose = true;
        } else if (token.startsWith("maxnoprogress=")) {
          heatConfig.resetAfterNoProgress = true;
          heatConfig.resetNoProgressIters = Integer.parseInt(token.substring(14).trim());
        } else if (token.equalsIgnoreCase("resettorandom")) {
          heatConfig.resetToRandom = true;
        } else {
          throw new RuntimeException("Unknown option token " + token + " for local heat pulses.");
        }
      }

      return new LocalHeatLocOpt(config, n, heatConfig);
    }

    final Newton newton = mapStringToLocOptMethod(input, config);
    if (newton == null) {
      return null;
    }

    final GenericLocOpt<Molecule, Geometry> fit = new NewtonAdaptor(newton);

    return fit;
  }

  @Override
  public GenericBackend<Molecule, Geometry> parseBackend(final String backend) throws Exception {

    final BackendFactory backFactory = new BackendFactory(config, params);

    return backFactory.parseBackend(backend);
  }

  public static Newton translateLocalOpt(final GenericLocOpt<Molecule, Geometry> locopt)
      throws Exception {

    if (locopt instanceof NewtonAdaptor) {
      // extract the newton part and return it
      final NewtonAdaptor adap = (NewtonAdaptor) locopt;
      return adap.getNewton();
    } else {
      return new InverseNewtonAdapter(locopt);
    }
  }

  /**
   * Maps the local optimization method given as a string into an integer array of optimization
   * methods.
   *
   * @param locMethod the string as found in the input describing the locopt choice
   * @param config the global configuration
   * @return Local optimization object.
   * @throws java.lang.Exception in case anything goes wrong. Important to catch and process (by
   *     printing and exiting).
   */
  public static Newton mapStringToLocOptMethod(final String locMethod, final GlobalConfig config)
      throws Exception {

    final String locOptString = locMethod.trim();
    Newton newton;
    if (locOptString.equalsIgnoreCase("external")) {
      newton = new ExternalLocOpt(config, 0, -FixedValues.NONCONVERGEDENERGY); // XXX better this!
    } else if (locOptString.startsWith("generic")) {
      newton = new GenericInterface(config);
    } else if (locOptString.startsWith("dummy")) {
      newton = new DummyLocOpt(config);
    } else if (locOptString.startsWith("adf:")) {
      final String inputStub = locOptString.substring(4).trim();
      newton = new ADFCaller(config, inputStub);
    } else if (locOptString.startsWith("xtb:")) {
      final String inputOpts = locOptString.substring(4).trim();

      XTBCaller.METHOD meth = XTBCaller.METHOD.GFN2XTB;
      XTBCaller.OPTLEVEL opt = XTBCaller.OPTLEVEL.NORMAL;
      boolean setEnvironment = true;
      String xControlFileName = null;
      boolean forceDelete = false;

      final String[] spl = inputOpts.split("\\,");
      for (final String s : spl) {
        if (s.trim().startsWith("method=")) {
          final String m = s.trim().substring("method=".length()).trim();
          if (m.equalsIgnoreCase("gfn-ff")) {
            meth = XTBCaller.METHOD.GFNFF;
          } else if (m.equalsIgnoreCase("gfn0-xtb")) {
            meth = XTBCaller.METHOD.GFN0XTB;
          } else if (m.equalsIgnoreCase("gfn1-xtb")) {
            meth = XTBCaller.METHOD.GFN1XTB;
          } else if (m.equalsIgnoreCase("gfn2-xtb")) {
            meth = XTBCaller.METHOD.GFN2XTB;
          } else {
            throw new RuntimeException("Illegal method " + m + " for XTB caller.");
          }
        } else if (s.trim().startsWith("optlevel=")) {
          final String o = s.trim().substring("optlevel=".length()).trim();
          if (o.equalsIgnoreCase("crude")) {
            opt = XTBCaller.OPTLEVEL.CRUDE;
          } else if (o.equalsIgnoreCase("sloppy")) {
            opt = XTBCaller.OPTLEVEL.SLOPPY;
          } else if (o.equalsIgnoreCase("loose")) {
            opt = XTBCaller.OPTLEVEL.LOOSE;
          } else if (o.equalsIgnoreCase("normal")) {
            opt = XTBCaller.OPTLEVEL.NORMAL;
          } else if (o.equalsIgnoreCase("tight")) {
            opt = XTBCaller.OPTLEVEL.TIGHT;
          } else if (o.equalsIgnoreCase("verytight")) {
            opt = XTBCaller.OPTLEVEL.VERYTIGHT;
          } else {
            throw new RuntimeException("Illegal optimization level " + o + " for XTB caller");
          }
        } else if (s.trim().equalsIgnoreCase("dontserialize")) {
          setEnvironment = false;
        } else if (s.trim().startsWith("xcontrol=")) {
          final String x = s.trim().substring("xcontrol=".length()).trim();
          xControlFileName = x;
          final File xControlFile = new File(xControlFileName);
          if (!xControlFile.isFile() || !xControlFile.canRead()) {
            throw new RuntimeException(
                "xcontrol file " + xControlFileName + " not found by XTB caller");
          }
        } else if (s.trim().equalsIgnoreCase("forcedelete")) {
          forceDelete = true;
        } else {
          throw new RuntimeException("Illegal XTB option " + s);
        }
      }

      newton = new XTBCaller(config, meth, opt, setEnvironment, xControlFileName, forceDelete);
    } else if (locOptString.startsWith("molpro:")) {
      String sTemp3 = locOptString.substring(7).trim();
      if (sTemp3.equalsIgnoreCase("am1/vdz") || sTemp3.equalsIgnoreCase("am1")) {
        newton = new MolproCaller(config, MolproCaller.METHOD.AM1);
      } else if (sTemp3.equalsIgnoreCase("hf/vdz")) {
        newton = new MolproCaller(config, MolproCaller.METHOD.HFVDZ);
      } else if (sTemp3.equalsIgnoreCase("b86/vdz")) {
        newton = new MolproCaller(config, MolproCaller.METHOD.B86VDZ);
      } else if (sTemp3.equalsIgnoreCase("mp2/avtz")) {
        newton = new MolproCaller(config, MolproCaller.METHOD.MP2AVTZ);
      } else if (sTemp3.equalsIgnoreCase("custom")) {
        newton = new MolproCaller(config, MolproCaller.METHOD.CUSTOM);
      } else if (sTemp3.equalsIgnoreCase("custom,nolocopt")) {
        newton = new MolproCaller(config, MolproCaller.METHOD.CUSTOMNOLOCOPT);
      } else {
        System.err.println("Wrong input to configure MOLPRO: " + sTemp3 + " using b86/vdz.");
        newton = new MolproCaller(config, MolproCaller.METHOD.B86VDZ);
      }
    } else if (locOptString.startsWith("mopac:")) {
      String sTemp3 = locOptString.substring(6).trim();
      if (sTemp3.equalsIgnoreCase("custom")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.CUSTOM);
      } else if (sTemp3.equalsIgnoreCase("mndo")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.MNDO);
      } else if (sTemp3.equalsIgnoreCase("am1")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.AM1);
      } else if (sTemp3.equalsIgnoreCase("pm3")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM3);
      } else if (sTemp3.equalsIgnoreCase("pm5")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM5);
      } else if (sTemp3.equalsIgnoreCase("pm6")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM6);
      } else if (sTemp3.equalsIgnoreCase("rm1")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.RM1);
      } else if (sTemp3.equalsIgnoreCase("pm7")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM7);
      } else if (sTemp3.equalsIgnoreCase("mndo,cosmo(ethanol)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.MNDO, 2.85f, 24.55f);
      } else if (sTemp3.equalsIgnoreCase("am1,cosmo(ethanol)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.AM1, 2.85f, 24.55f);
      } else if (sTemp3.equalsIgnoreCase("pm3,cosmo(ethanol)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM3, 2.85f, 24.55f);
      } else if (sTemp3.equalsIgnoreCase("pm5,cosmo(ethanol)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM5, 2.85f, 24.55f);
      } else if (sTemp3.equalsIgnoreCase("pm6,cosmo(ethanol)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM6, 2.85f, 24.55f);
      } else if (sTemp3.equalsIgnoreCase("rm1,cosmo(ethanol)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.RM1, 2.85f, 24.55f);
      } else if (sTemp3.equalsIgnoreCase("pm7,cosmo(ethanol)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM7, 2.85f, 24.55f);
      } else if (sTemp3.equalsIgnoreCase("mndo,cosmo(water)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.MNDO, 1.93f, 78.39f);
      } else if (sTemp3.equalsIgnoreCase("am1,cosmo(water)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.AM1, 1.93f, 78.39f);
      } else if (sTemp3.equalsIgnoreCase("pm3,cosmo(water)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM3, 1.93f, 78.39f);
      } else if (sTemp3.equalsIgnoreCase("pm5,cosmo(water)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM5, 1.93f, 78.39f);
      } else if (sTemp3.equalsIgnoreCase("pm6,cosmo(water)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM6, 1.93f, 78.39f);
      } else if (sTemp3.equalsIgnoreCase("rm1,cosmo(water)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.RM1, 1.93f, 78.39f);
      } else if (sTemp3.equalsIgnoreCase("pm7,cosmo(water)")) {
        newton = new MopacCaller(config, MopacCaller.METHOD.PM7, 1.93f, 78.39f);
      } else {
        System.err.println("Wrong input to configure MOPAC: " + sTemp3 + " using pm3.");
        newton = new MopacCaller(config, MopacCaller.METHOD.PM3);
      }
    } else if (locOptString.startsWith("dftb+:")) {
      String sTemp3 = locOptString.substring(6).trim();
      if (sTemp3.startsWith("steepest")) {
        String sTemp4 = sTemp3.substring(8).trim();
        if (sTemp4.equalsIgnoreCase(",scc")) {
          newton =
              new DFTBplusCaller(config, DFTBplusCaller.WHICHDRIVER.steepestdescent, true, false);
        } else {
          newton =
              new DFTBplusCaller(config, DFTBplusCaller.WHICHDRIVER.steepestdescent, false, false);
        }
      } else if (sTemp3.startsWith("conjugate")) {
        String sTemp4 = sTemp3.substring(9).trim();
        if (sTemp4.equalsIgnoreCase(",scc")) {
          newton =
              new DFTBplusCaller(config, DFTBplusCaller.WHICHDRIVER.conjugategradient, true, false);
        } else {
          newton =
              new DFTBplusCaller(
                  config, DFTBplusCaller.WHICHDRIVER.conjugategradient, false, false);
        }
      } else if (sTemp3.startsWith("custom")) {
        newton =
            new DFTBplusCaller(config, DFTBplusCaller.WHICHDRIVER.conjugategradient, true, true);
      } else {
        System.err.println(
            "Wrong input to configure DFTB+: " + sTemp3 + " using conjugate gradient with SCC.");
        newton =
            new DFTBplusCaller(config, DFTBplusCaller.WHICHDRIVER.conjugategradient, true, false);
      }
    } else if (locOptString.startsWith("orca:")) {
      String sTemp3 = locOptString.substring(5).trim();
      final String auxFile = config.outputFolder + "-orca.aux";
      if (sTemp3.equalsIgnoreCase("mndo")) {
        newton = new OrcaCaller(config, 0, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("am1")) {
        newton = new OrcaCaller(config, 1, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("pm3")) {
        newton = new OrcaCaller(config, 2, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("b3lyp/vdz")) {
        newton = new OrcaCaller(config, 3, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("bp86/svp")) {
        newton = new OrcaCaller(config, 4, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("bp86/tzvp")) {
        newton = new OrcaCaller(config, 5, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("ubp86/tzvp")) {
        newton = new OrcaCaller(config, 6, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("b3lyp/6-31+g**")) {
        newton = new OrcaCaller(config, 7, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("b3lyp/6-31+g**,vdw")) {
        newton = new OrcaCaller(config, 8, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("mp2/aug-cc-pVDZ")) {
        newton = new OrcaCaller(config, 40, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("ccsd(t)/tzvpp")) {
        newton = new OrcaCaller(config, 50, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("custom")) {
        newton = new OrcaCaller(config, -1, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("custom,nolocopt")) {
        newton = new OrcaCaller(config, -2, false, auxFile);
      } else if (sTemp3.equalsIgnoreCase("mndo,fileout")) {
        newton = new OrcaCaller(config, 0, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("am1,fileout")) {
        newton = new OrcaCaller(config, 1, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("pm3,fileout")) {
        newton = new OrcaCaller(config, 2, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("b3lyp/vdz,fileout")) {
        newton = new OrcaCaller(config, 3, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("bp86/svp,fileout")) {
        newton = new OrcaCaller(config, 4, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("bp86/tzvp,fileout")) {
        newton = new OrcaCaller(config, 5, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("ubp86/tzvp,fileout")) {
        newton = new OrcaCaller(config, 6, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("b3lyp/6-31+g**,fileout")) {
        newton = new OrcaCaller(config, 7, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("b3lyp/6-31+g**,vdw,fileout")) {
        newton = new OrcaCaller(config, 8, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("mp2/aug-cc-pVDZ,fileout")) {
        newton = new OrcaCaller(config, 40, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("ccsd(t)/tzvpp,fileout")) {
        newton = new OrcaCaller(config, 50, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("custom,fileout")) {
        newton = new OrcaCaller(config, -1, true, auxFile);
      } else if (sTemp3.equalsIgnoreCase("custom,nolocopt,fileout")) {
        newton = new OrcaCaller(config, -2, true, auxFile);
      } else {
        throw new RuntimeException("Wrong input to configure ORCA: " + sTemp3 + ".");
      }
    } else if (locOptString.startsWith("tinker:")) {
      String sTemp3 = locOptString.substring(7).trim();
      if (sTemp3.equalsIgnoreCase("minimize")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.minimize, false, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("newton")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.newton, false, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("optimize")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.optimize, false, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("custom")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("minimize,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.minimize, true, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("newton,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.newton, true, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("optimize,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.optimize, true, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("custom,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, true, false, false, false);
      } else if (sTemp3.equalsIgnoreCase("minimize,fileout")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.minimize, false, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("newton,fileout")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.newton, false, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("optimize,fileout")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.optimize, false, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("custom,fileout")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("minimize,fileout,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.minimize, true, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("newton,fileout,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.newton, true, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("optimize,fileout,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.optimize, true, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("custom,fileout,solvate")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, true, true, false, false);
      } else if (sTemp3.equalsIgnoreCase("ubercustom")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, false, true, false);
      } else if (sTemp3.equalsIgnoreCase("ubercustom,fileout")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, true, true, false);
      } else if (sTemp3.equalsIgnoreCase("ubercustom2")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, false, true, true);
      } else if (sTemp3.equalsIgnoreCase("ubercustom2,fileout")) {
        newton = new TinkerCaller(config, TinkerCaller.METHOD.custom, false, true, true, true);
      } else {
        System.err.println("Wrong input to configure tinker: " + sTemp3 + " using minimize.");
        newton = new TinkerCaller(config, TinkerCaller.METHOD.minimize, false, false, false, false);
      }
    } else if (locOptString.startsWith("tinkermd:")) {
      String sTemp3 = locOptString.substring(9).trim();
      if (sTemp3.equalsIgnoreCase("pss")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.pss, false);
      } else if (sTemp3.equalsIgnoreCase("anneal")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.anneal, false);
      } else if (sTemp3.equalsIgnoreCase("monte")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.monte, false);
      } else if (sTemp3.equalsIgnoreCase("sniffer")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.sniffer, false);
      } else if (sTemp3.equalsIgnoreCase("pss,solvate")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.pss, true);
      } else if (sTemp3.equalsIgnoreCase("anneal,solvate")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.anneal, true);
      } else if (sTemp3.equalsIgnoreCase("monte,solvate")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.monte, true);
      } else if (sTemp3.equalsIgnoreCase("sniffer,solvate")) {
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.sniffer, true);
      } else {
        System.err.println("Wrong input to configure tinkermd: " + sTemp3 + " using pss.");
        newton = new TinkerMDCaller(config, TinkerMDCaller.METHOD.pss, false);
      }
    } else if (locOptString.startsWith("gulp:")) {
      final String ankString = locOptString.substring(5);
      newton = new GulpCaller(config, "gulp-aux.aux", ankString);
    } else if (locOptString.startsWith("openbabel:")) {
      String sTemp3 = locOptString.substring(10).trim();
      if (sTemp3.equalsIgnoreCase("ghemical")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.GHEMICAL, false);
      } else if (sTemp3.equalsIgnoreCase("mmff94")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.MMFF94, false);
      } else if (sTemp3.equalsIgnoreCase("mmff94s")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.MMFF94s, false);
      } else if (sTemp3.equalsIgnoreCase("uff")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.UFF, false);
      } else if (sTemp3.equalsIgnoreCase("ghemical,nocache")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.GHEMICAL, true);
      } else if (sTemp3.equalsIgnoreCase("mmff94,nocache")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.MMFF94, true);
      } else if (sTemp3.equalsIgnoreCase("mmff94s,nocache")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.MMFF94s, true);
      } else if (sTemp3.equalsIgnoreCase("uff,nocache")) {
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.UFF, true);
      } else {
        System.err.println(
            "Wrong input to configure OpenBabel: "
                + sTemp3
                + " using the GHEMICAL force field w/ caching.");
        newton = new OpenBabelCaller(config, OpenBabelCaller.FORCEFIELD.GHEMICAL, false);
      }
    } else if (locOptString.startsWith("amber")) {
      String sTemp3 = locOptString.substring(5).trim();
      if (sTemp3.equalsIgnoreCase("")) {
        newton = new AmberCaller(config, false);
      } else if (sTemp3.equalsIgnoreCase(":implicitwater")) {
        newton = new AmberCaller(config, true);
      } else {
        System.err.println(
            "Wrong input to configure amber: " + sTemp3 + " using no implicit water.");
        newton = new AmberCaller(config, false);
      }
    } else if (locOptString.startsWith("gromacs")) {
      final String sTemp3 = locOptString.substring(7).trim();
      if (sTemp3.equalsIgnoreCase("")) {
        newton = new GromacsCaller(config, GromacsCaller.QMENGINE.NONE);
      } else if (sTemp3.equalsIgnoreCase(":mopac")) {
        newton = new GromacsCaller(config, GromacsCaller.QMENGINE.MOPAC);
      } else if (sTemp3.equalsIgnoreCase(":orca")) {
        newton = new GromacsCaller(config, GromacsCaller.QMENGINE.Orca);
      } else if (sTemp3.equalsIgnoreCase(":gaussian")) {
        newton = new GromacsCaller(config, GromacsCaller.QMENGINE.Gaussian);
      } else if (sTemp3.equalsIgnoreCase(":gamess-uk")) {
        newton = new GromacsCaller(config, GromacsCaller.QMENGINE.GamessUK);
      } else {
        System.err.println("Wrong input to configure gromacs: " + sTemp3 + " using pure MM now.");
        newton = new GromacsCaller(config, GromacsCaller.QMENGINE.NONE);
      }
    } else if (locOptString.startsWith("namd")) {
      String sTemp3 = locOptString.substring(4).trim();
      if (sTemp3.equalsIgnoreCase("")) {
        newton = new NAMDCaller(config);
      } else {
        System.err.println(
            "Wrong input to configure namd: " + sTemp3 + " using no implicit water.");
        newton = new NAMDCaller(config);
      }
    } else if (locOptString.startsWith("penaltylocopt:")) {
      String sTemp3 = locOptString.substring(14).trim();

      final Newton newt1 = mapStringToLocOptMethod(sTemp3, config);
      newton = new PenaltyLocOpt(config, newt1);
    } else if (locOptString.startsWith("mndo:")) {
      String sTemp3 = locOptString.substring(5).trim();

      if (sTemp3.equalsIgnoreCase("mndod")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MNDOd, false);
      } else if (sTemp3.equalsIgnoreCase("om3")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.OM3, false);
      } else if (sTemp3.equalsIgnoreCase("pm3")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.PM3, false);
      } else if (sTemp3.equalsIgnoreCase("om2")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.OM2, false);
      } else if (sTemp3.equalsIgnoreCase("om1")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.OM1, false);
      } else if (sTemp3.equalsIgnoreCase("am1")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.AM1, false);
      } else if (sTemp3.equalsIgnoreCase("mndoc")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MNDOC, false);
      } else if (sTemp3.equalsIgnoreCase("mndo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MNDO, false);
      } else if (sTemp3.equalsIgnoreCase("mindo3")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MINDO3, false);
      } else if (sTemp3.equalsIgnoreCase("cndo2")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.CNDO2, false);
      } else if (sTemp3.equalsIgnoreCase("scc-dftb")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.SCCDFTB, false);
      } else if (sTemp3.equalsIgnoreCase("scc-dftbj")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.SCCDFTBJorgensen, false);
      } else if (sTemp3.equalsIgnoreCase("mndod:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MNDOd, true);
      } else if (sTemp3.equalsIgnoreCase("om3:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.OM3, true);
      } else if (sTemp3.equalsIgnoreCase("pm3:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.PM3, true);
      } else if (sTemp3.equalsIgnoreCase("om2:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.OM2, true);
      } else if (sTemp3.equalsIgnoreCase("om1:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.OM1, true);
      } else if (sTemp3.equalsIgnoreCase("am1:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.AM1, true);
      } else if (sTemp3.equalsIgnoreCase("mndoc:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MNDOC, true);
      } else if (sTemp3.equalsIgnoreCase("mndo:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MNDO, true);
      } else if (sTemp3.equalsIgnoreCase("mindo3:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.MINDO3, true);
      } else if (sTemp3.equalsIgnoreCase("cndo2:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.CNDO2, true);
      } else if (sTemp3.equalsIgnoreCase("scc-dftb:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.SCCDFTB, true);
      } else if (sTemp3.equalsIgnoreCase("scc-dftbj:cosmo")) {
        newton = new MNDOCaller(config, MNDOCaller.METHOD.SCCDFTBJorgensen, true);
      } else {
        System.err.println("Wrong input to configure MNDO: " + sTemp3 + " using OM3.");
        newton = new MNDOCaller(config, MNDOCaller.METHOD.OM3, false);
      }

    } else if (locOptString.startsWith("cp2k:")) {
      final String sTemp = locOptString.substring(5).trim();
      if (sTemp.equalsIgnoreCase("B3LYP")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.B3LYP);
      } else if (sTemp.equalsIgnoreCase("BLYP")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.BLYP);
      } else if (sTemp.equalsIgnoreCase("BP")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.BP);
      } else if (sTemp.equalsIgnoreCase("HCTH120")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.HCTH120);
      } else if (sTemp.equalsIgnoreCase("NONE")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.NONE);
      } else if (sTemp.equalsIgnoreCase("NO_SHORTCUT")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.NO_SHORTCUT);
      } else if (sTemp.equalsIgnoreCase("OLYP")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.OLYP);
      } else if (sTemp.equalsIgnoreCase("PADE")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.PADE);
      } else if (sTemp.equalsIgnoreCase("PBE")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.PBE);
      } else if (sTemp.equalsIgnoreCase("PBE0")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.PBE0);
      } else if (sTemp.equalsIgnoreCase("TPSS")) {
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.TPSS);
      } else {
        System.err.println("Wrong input to configure CP2K: " + sTemp + " using BLYP.");
        newton = new CP2KCaller(config, CP2KCaller.FUNCTIONAL.BLYP);
      }
    } else if (locOptString.startsWith("cpmd:")) {
      final String sTemp = locOptString.substring(5).trim();
      if (sTemp.equalsIgnoreCase("SONLY")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.SONLY);
      } else if (sTemp.equalsIgnoreCase("LDA")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.LDA);
      } else if (sTemp.equalsIgnoreCase("BONLY")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.BONLY);
      } else if (sTemp.equalsIgnoreCase("BP")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.BP);
      } else if (sTemp.equalsIgnoreCase("BLYP")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.BLYP);
      } else if (sTemp.equalsIgnoreCase("XLYP")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.XLYP);
      } else if (sTemp.equalsIgnoreCase("GGA")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.GGA);
      } else if (sTemp.equalsIgnoreCase("PBE")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.PBE);
      } else if (sTemp.equalsIgnoreCase("PBES")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.PBES);
      } else if (sTemp.equalsIgnoreCase("REVPBE")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.REVPBE);
      } else if (sTemp.equalsIgnoreCase("HCTH")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.HCTH);
      } else if (sTemp.equalsIgnoreCase("OPTX")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.OPTX);
      } else if (sTemp.equalsIgnoreCase("OLYP")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.OLYP);
      } else if (sTemp.equalsIgnoreCase("TPSS")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.TPSS);
      } else if (sTemp.equalsIgnoreCase("PBE0")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.PBE0);
      } else if (sTemp.equalsIgnoreCase("B1LYP")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.B1LYP);
      } else if (sTemp.equalsIgnoreCase("B3LYP")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.B3LYP);
      } else if (sTemp.equalsIgnoreCase("X3LYP")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.X3LYP);
      } else if (sTemp.equalsIgnoreCase("CUSTOM")) {
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.CUSTOM);
      } else {
        System.err.println("Wrong input to configure CP2K: " + sTemp + " using BLYP.");
        newton = new CPMDCaller(config, CPMDCaller.FUNCTIONAL.BLYP);
      }
    } else if (locOptString.startsWith("lotriff")) {
      newton = new LotriffCaller(config);
    } else if (locOptString.startsWith("vasp:")) {
      final String s = locOptString.substring(5).trim();
      final String[] sa = s.split(";");
      final boolean symmetrize = (sa.length >= 3) ? Boolean.parseBoolean(sa[2]) : false;
      final double height =
          (sa.length >= 4) ? Double.parseDouble(sa[3]) * Constants.ANGTOBOHR : 0.0;
      int whichref = -1;
      String refAtom = "";
      double[] refCoords = null;
      if (sa.length >= 5) {
        final String[] sa2 = sa[4].trim().split("\\,");
        whichref = Integer.parseInt(sa2[0].trim());
        refAtom = sa2[1].trim();
        refCoords = new double[3];
        for (int c = 0; c < 3; c++) {
          refCoords[c] = Double.parseDouble(sa2[2 + c].trim()) * Constants.ANGTOBOHR;
        }
      }
      newton =
          new VASPCaller(config, sa[0], sa[1], symmetrize, height, whichref, refAtom, refCoords);
    } else {
      return null;
    }

    return newton;
  }
}
