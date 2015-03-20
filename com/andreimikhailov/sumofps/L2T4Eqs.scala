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

/** This is to show that there is no cohomology of the type lambda2theta4
  *
  */
object L2T4Eqs extends App {

  println("=== calculating LMT4 ===")
  val LMT4a = capless(
    unknown(1) *: word("l^m_t^t_t^t") +
      unknown(2) *: word("l^t_m^t_t^t") +
      unknown(3) *: word("l^t_t^m_t^t") +
      unknown(4) *: word("l^t_t^t_m^t") +
      unknown(5) *: word("l^t_t^t_t^m") +
      unknown(6) *: word("t^l_m^t_t^t") +
      unknown(7) *: word("t^l_t^m_t^t") +
      unknown(8) *: word("t^l_t^t_m^t") +
      unknown(9) *: word("t^l_t^t_t^m") +
      unknown(10) *: word("t^t_l^m_t^t") +
      unknown(11) *: word("t^t_l^t_m^t") +
      unknown(12) *: word("t^t_l^t_t^m") +
      unknown(13) *: word("t^t_t^l_m^t") +
      unknown(14) *: word("t^t_t^l_t^m") +
      unknown(15) *: word("t^t_t^t_l^m") +
      unknown(16) *: ((word("t^t_t^t_").trace) *: word("l^m")) +
      unknown(17) *: ((word("l^t_m^t_").trace) *: word("t^t")) +
      unknown(18) *: ((word("l^t_t^t_").trace) *: word("m^t")) +
      unknown(19) *: ((word("l^t_t^t_").trace) *: word("t^m")) +
      (unknown(20) * (word("l^t_").trace) * (word("m^t_").trace)) *: word("t^t") +
      (unknown(21) * (word("l^t_").trace)) *: word("m^t_t^t") +
      (unknown(22) * (word("l^t_").trace)) *: word("t^m_t^t") +
      (unknown(23) * (word("l^t_").trace)) *: word("t^t_m^t") +
      (unknown(24) * (word("l^t_").trace)) *: word("t^t_t^m")
  )

  val LMT4b = capless(
    unknown(1) *: word("m^l_t^t_t^t") +
      unknown(2) *: word("m^t_l^t_t^t") +
      unknown(3) *: word("m^t_t^l_t^t") +
      unknown(4) *: word("m^t_t^t_l^t") +
      unknown(5) *: word("m^t_t^t_t^l") +
      unknown(6) *: word("t^m_l^t_t^t") +
      unknown(7) *: word("t^m_t^l_t^t") +
      unknown(8) *: word("t^m_t^t_l^t") +
      unknown(9) *: word("t^m_t^t_t^l") +
      unknown(10) *: word("t^t_m^l_t^t") +
      unknown(11) *: word("t^t_m^t_l^t") +
      unknown(12) *: word("t^t_m^t_t^l") +
      unknown(13) *: word("t^t_t^m_l^t") +
      unknown(14) *: word("t^t_t^m_t^l") +
      unknown(15) *: word("t^t_t^t_m^l") +
      unknown(16) *: ((word("t^t_t^t_").trace) *: word("m^l")) +
      unknown(17) *: ((word("m^t_l^t_").trace) *: word("t^t")) +
      unknown(18) *: ((word("m^t_t^t_").trace) *: word("l^t")) +
      unknown(19) *: ((word("m^t_t^t_").trace) *: word("t^l")) +
      (unknown(20) * (word("m^t_").trace) * (word("l^t_").trace)) *: word("t^t") +
      (unknown(21) * (word("m^t_").trace)) *: word("l^t_t^t") +
      (unknown(22) * (word("m^t_").trace)) *: word("t^l_t^t") +
      (unknown(23) * (word("m^t_").trace)) *: word("t^t_l^t") +
      (unknown(24) * (word("m^t_").trace)) *: word("t^t_t^l")
  )

  println("=== calculating Q(LMT4a + LMT4b) ===")

  /** Calculating the coefficients in the ansatz of the form
    * <br>
    * <img src="http://andreimikhailov.com/media/minitheta-scaladoc/TLTTLTplusTMTTMT.png">
    */

  val zmapSymm =
    Equations(
      collectEquationsFromMat(
        QLmat(capless(word("t^m_t^t_m^t"))) +
          QRmat(LMT4a + LMT4b) + QLmat(LMT4a + LMT4b) +
          QRmat(capless(word("t^l_t^t_l^t")))),
      varList(0, 300)
    ).solve

  println("=== calculating Q(LMT4a - LMT4b) ===")
  /** Calculating the coefficients in the ansatz of the form
    * <br>
    * <img src="http://andreimikhailov.com/media/minitheta-scaladoc/TLTTLTminusTMTTMT.png">
    */

  val zmapAntiSymm =
    Equations(
      collectEquationsFromMat(
        QLmat(capless(word("t^m_t^t_m^t"))) +
          QRmat(LMT4a + (scalarESCA(-1) *: LMT4b)) +
          QLmat(LMT4a + (scalarESCA(-1) *: LMT4b)) +
          QRmat(scalarESCA(-1) *: capless(word("t^l_t^t_l^t")))),
      varList(0, 300)
    ).solve

  println("Symmetric solutions: " + zmapSymm)
  println("Antisymmetric solutions: " + zmapAntiSymm)
  println("Since both are None, we conclude that there is no cohomology of this type.")
}