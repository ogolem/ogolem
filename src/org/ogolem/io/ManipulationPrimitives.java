/**
Copyright (c) 2012-2014, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
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
package org.ogolem.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import org.ogolem.helpers.Tuple3D;

/**
 * A set of manipulation primitives.
 * @author Johannes Dieterich
 * @version 2020-04-29
 */
public class ManipulationPrimitives {
	
	private static final boolean DEBUG = false;
    
    public static void remove(final String path) throws IOException{
        
        final File f = new File(path);
        // if this is a file, no recursive call is needed
        if(f.isFile()){
            if(!f.delete()){ System.err.println("ERROR: Couldn't delete file " + path);}
            return;
        }
        
        // handle the directory case
        final File[] files = f.listFiles();
        for (final File file : files) {
            if (file.isDirectory()) {
                // recursive call... :-)
                remove(file.getAbsolutePath());
            } else {
                if (!file.delete()) {
                    System.err.println("ERROR: File couldn't be deleted. " + file);
                }
            }
        }
        if (!f.delete()) {
            // no success. for whatever reason.
            throw new IOException("No success in deleting directory " + path + " .");
        }
    }
    
    public static void moveFile(final String original, final String moved) throws IOException {
        final File f = new File(original);
        final boolean sucess = f.renameTo(new File(moved));
        if (!sucess) {
            throw new IOException("Couldn't move file " + original + " to " + moved + ".");
        }
    }
    
    /**
     * Moves old folders up to a certain maximum to keep that data.
     * @param folderName The folder name that is supposed to be moved.
     * @param maxFolders the number of folders we want to store
     * @param pathPrefix the full path prefix. e.g. user.dir
     * @throws java.io.IOException if something does wrong with the I/O
     */
    public static void moveOldFolders(final String folderName, final int maxFolders,
            final String pathPrefix) throws IOException {
            	
        System.out.println("Working on " + folderName);
                
        String fullPath = pathPrefix;
        if (!fullPath.endsWith(System.getProperty("file.separator"))) {
            fullPath += System.getProperty("file.separator");
        }
        
        if(folderName.startsWith(pathPrefix)){
            // glue together, folder name is relative path
            fullPath += folderName;
        } else {
            // just use the folder name, absolute path
            fullPath = folderName;
        }

        int i = maxFolders;
        do {
        	
        	if(DEBUG) System.out.println("Moving for #" + i  + " and " + fullPath);
        	
            if (i == maxFolders) {
                final String sFolderPath = fullPath + "." + i;
                final File f = new File(sFolderPath);
                if (f.exists()) {
                	if(DEBUG) System.out.println("Removing " + sFolderPath);
                    // if folder.$iMax (here: maxFolders) exists, then delete it and its content.
                    remove(sFolderPath);
                }

            } else if (i == 0) {
                final File f = new File(fullPath);
                if (f.exists()) {
                    // move folder away to folder.1
                    final String sFolderNew = fullPath + "." + (i + 1);
                    final File fNew = new File(sFolderNew);
                    if(DEBUG) System.out.println("Moving " + fullPath + " to " + sFolderNew);
                    if (!f.renameTo(fNew)) {
                        // not successful. can have various reasons according to the documentation.
                        // for example an existing folder of this name makes the moving impossible
                        // but that is handled with the other else if.
                        throw new IOException("Moving the old folder from " + fullPath + " to " + sFolderNew + " failed.");
                    }
                }
            } else {
                final String oldFolderPath = fullPath + "." + i;
                final File f = new File(oldFolderPath);
                // for folders i = 1 to max (exclusive)
                if (f.exists()) {
                    // move folder to the one higher number
                    final String folderPathNew = fullPath + "." + (i + 1);
                    final File fNew = new File(folderPathNew);
                    if(DEBUG) System.out.println("Moving " + fullPath + " to " + oldFolderPath);
                    if (!f.renameTo(fNew)) {
                        // not successful.
                        throw new IOException("Moving the old folder from " + oldFolderPath + " to " + folderPathNew + " failed.");
                    }
                }
            }

            i--;
        } while (i >= 0);
    }
    
    /**
     * Marks a file executable.
     * @param fileName the file of the file supposed to be marked executable
     * @throws Exception if something goes wrong
     */
    public static final void markExecutable(final String fileName) throws Exception {
        
        final File f = new File(fileName);
        f.setExecutable(true);
    }
    
    /**
     * Obtains the input directory, output directory, and base name for outputs from the input file name.
     * @param inputFile
     * @return a Tuple containing the input directory (absolute path), output directory (relative path), and base name (relative path) in that order.
     */
    public static Tuple3D<String, String,String> outDirAndBaseName(final String inputFile){
        
    	final File f = new File(inputFile);
    	final String absolutePath = f.getAbsolutePath();
    	
    	if(DEBUG) System.out.println("Absolute path" + absolutePath);
    	
    	final String parentFolder = Paths.get(absolutePath).getParent().toString();
        final String fileName = absolutePath;
        final int indexOfSuffix = fileName.indexOf(".ogo");
        
        final String basePath = fileName.substring(0,indexOfSuffix);
        final String split[] = basePath.split("\\" + File.separator);
        final String baseName = split[split.length-1]; // last entry is the real base
        final String outFolder = baseName;//parentFolder + File.separator + baseName;
        
        if(DEBUG) {
        	System.out.println("Returning following paths:");
        	System.out.println("parentFolder " + parentFolder);
        	System.out.println("outFolder " + outFolder);
        	System.out.println("baseName " + baseName);
        }
        
        return new Tuple3D<>(parentFolder, outFolder, baseName);
    }
}
