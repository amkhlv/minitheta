/*
 * minitheta
 * a Scala library for doing symbolic computations with polynomials of Grassmann variables
 *
 * Author: Andrei Mikhailov <a.mkhlv@gmail.com>
 *
 * Copyright (C) 2012
 *
 * This program is free software which I release under the GNU General Public
 * License. You may redistribute and/or modify this program under the terms
 * of that license as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * To get a copy of the GNU General Puplic License,  write to the
 * Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

package com.andreimikhailov.sumofps

import com.andreimikhailov.minitheta.core._
import com.andreimikhailov.brst._

object L2T2Eqs extends App {
  println("=== calculating LMT4 ===")

  /** This is the main ansatz for the expression of the form $\lambda_L\lambda_R\theta\theta$
    * Notice that it is parametrized by 6 + 6 + 5 = 17 coefficients
    */
  val matr = capless(
    matScalarUnknown(11) * word("l^m_t^t") +
      matScalarUnknown(12) * word("l^t_m^t") +
      matScalarUnknown(13) * word("l^t_t^m") +
      matScalarUnknown(14) * word("t^l_m^t") +
      matScalarUnknown(15) * word("t^l_t^m") +
      matScalarUnknown(16) * word("t^t_l^m") +
      //
      matScalarUnknown(21) * word("m^l_t^t") +
      matScalarUnknown(22) * word("m^t_l^t") +
      matScalarUnknown(23) * word("m^t_t^l") +
      matScalarUnknown(24) * word("t^m_l^t") +
      matScalarUnknown(25) * word("t^m_t^l") +
      matScalarUnknown(26) * word("t^t_m^l") +
      //
      matScalarUnknown(1) * matScalar(word("l^m_").trace) * word("t^t") +
      matScalarUnknown(2) * matScalar(word("l^t_").trace) * word("m^t") +
      matScalarUnknown(3) * matScalar(word("l^t_").trace) * word("t^m") +
      matScalarUnknown(4) * matScalar(word("m^t_").trace) * word("l^t") +
      matScalarUnknown(5) * matScalar(word("m^t_").trace) * word("t^l")
  )

  /** Some of these expressions may be $Q_LQ_R$ of $\theta\theta\theta\theta$,
    * here is $\theta\theta\theta\theta$
    */
  val bndOf = capless(matScalarUnknown(100) * word("t^t_t^t"))

  val setOfVariables = varList(0, 300)
  val zSolutionForBeingZero = Equations(
    collectEquationsFromMat(matr), setOfVariables
  ).solve match {
    case Some(zmap) => zmap
    case None => throw new RuntimeException("NO MATCH")
  }

  val zSolutionForBeingAnnihilatedByQ = Equations(
    collectEquationsFromMat(QLmat(matr) + QRmat(matr)),
    setOfVariables
  ).solve match {
    case Some(zmap) => zmap
    case None => throw new RuntimeException("NO MATCH")
  }

  println("Notice that there are 17 coefficients.")
  println("But some of these expressions are identically zero.")
  println("The condition for the expression being zero is:")
  for (k <- zSolutionForBeingZero.keys) println(k + "--->  " + zSolutionForBeingZero(k))
  println("Altogether " + zSolutionForBeingZero.keys.toList.length + " constraints.")
  println("")
  println("And the conditions for being annihilated by Q are:")
  for (k <- zSolutionForBeingAnnihilatedByQ.keys) println(k + "--->  " + zSolutionForBeingAnnihilatedByQ(k))
  println("Altogether " + zSolutionForBeingAnnihilatedByQ.keys.toList.length +
    " constraints for the expression to be in ker(Q)")
  println("")
  println("And remember that 1 expression is $Q_LQ_R$ of $\\theta\\theta\\theta\\theta$")

}