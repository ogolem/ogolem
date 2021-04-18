/**
 * Copyright (c) 2009-2010, J. M. Dieterich and B. Hartke 2010, J. M. Dieterich 2020, J. M.
 * Dieterich and B. Hartke All rights reserved.
 *
 * <p>Redistribution and use in source and binary forms, with or without modification, are permitted
 * provided that the following conditions are met:
 *
 * <p>Redistributions of source code must retain the above copyright notice, this list of conditions
 * and the following disclaimer.
 *
 * <p>Redistributions in binary form must reproduce the above copyright notice, this list of
 * conditions and the following disclaimer in the documentation and/or other materials provided with
 * the distribution.
 *
 * <p>All advertising materials mentioning features or use of this software must display the
 * following acknowledgement:
 *
 * <p>This product includes software of the ogolem.org project developed by J. M. Dieterich and B.
 * Hartke (Christian-Albrechts-University Kiel, Germany) and contributors.
 *
 * <p>Neither the name of the ogolem.org project, the University of Kiel nor the names of its
 * contributors may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * <p>THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package org.ogolem.adaptive;

/**
 * An interface describing what an adaptive local optimization needs to know to work well both with
 * the core package and the adaptive scheme.
 *
 * @author Johannes Dieterich
 * @version 2020-12-29
 */
interface AdaptiveNewton extends org.ogolem.core.Newton, Adaptivable {

  @Override
  AdaptiveNewton copy();
}
