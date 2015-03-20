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


object L2T3Eqs extends App {

  println("=== calculating L2T3 ===")
  val L2T3 = unknown(201) *: (traceless(word("t^l_")) * word("t^l_t"))
  println("=== calculating M2T3 ===")
  val M2T3 = unknown(202) *: (traceless(word("t^m_")) * word("t^m_t"))

  println("=== running checks ===")
  val zeroLeft = QLmat(L2T3)
  val zeroRight = QRmat(M2T3)
  println(zeroLeft)
  println(zeroRight)
  assert(QLmat(L2T3) == mscalar(0, 1, 4))
  assert(QRmat(M2T3) == mscalar(0, 1, 4))

  println("=== calculating LMT3 ===")
  val LMT3 =
    unknown(1) *: word("l^m_t^t_t") +
      unknown(2) *: word("l^t_m^t_t") +
      unknown(3) *: word("l^t_t^m_t") +
      unknown(4) *: word("l^t_t^t_m") +
      unknown(5) *: word("t^l_m^t_t") +
      unknown(6) *: word("t^l_t^m_t") +
      unknown(7) *: word("t^l_t^t_m") +
      unknown(8) *: word("t^t_l^m_t") +
      unknown(9) *: word("t^t_l^t_m") +
      unknown(10) *: word("t^t_t^l_m") +
      unknown(11) *: ((word("l^m_").trace) *: word("t^t_t")) +
      unknown(12) *: ((word("l^t_").trace) *: word("m^t_t")) +
      unknown(13) *: ((word("l^t_").trace) *: word("t^m_t")) +
      unknown(14) *: ((word("l^t_").trace) *: word("t^t_m")) +
      (unknown(15) * (word("l^t_").trace) * (word("m^t_").trace)) *: word("t") +
      (unknown(16) * (word("l^m_t^t_").trace)) *: word("t") +
      (unknown(17) * (word("l^t_m^t_").trace)) *: word("t") +
      (unknown(18) * (word("l^t_t^m_").trace)) *: word("t") +
      (unknown(19) * (word("l^t_t^t_").trace)) *: word("m") +
      unknown(51) *: word("m^l_t^t_t") +
      unknown(52) *: word("m^t_l^t_t") +
      unknown(53) *: word("m^t_t^l_t") +
      unknown(54) *: word("m^t_t^t_l") +
      unknown(55) *: word("t^m_l^t_t") +
      unknown(56) *: word("t^m_t^l_t") +
      unknown(57) *: word("t^m_t^t_l") +
      unknown(58) *: word("t^t_m^l_t") +
      unknown(59) *: word("t^t_m^t_l") +
      unknown(60) *: word("t^t_t^m_l") +
      unknown(61) *: ((word("m^l_").trace) *: word("t^t_t")) +
      unknown(62) *: ((word("m^t_").trace) *: word("l^t_t")) +
      unknown(63) *: ((word("m^t_").trace) *: word("t^l_t")) +
      unknown(64) *: ((word("m^t_").trace) *: word("t^t_l")) +
      (unknown(65) * (word("m^t_").trace) * (word("l^t_").trace)) *: word("t") +
      (unknown(66) * (word("m^l_t^t_").trace)) *: word("t") +
      (unknown(67) * (word("m^t_l^t_").trace)) *: word("t") +
      (unknown(68) * (word("m^t_t^l_").trace)) *: word("t") +
      (unknown(69) * (word("m^t_t^t_").trace)) *: word("l")

  println("=== Calculating BRST transformation and extracting equations ===")

  val zmap =
    Equations(collectEquationsFromMat(
      QLmat(L2T3 + M2T3 + LMT3) + QRmat(L2T3 + M2T3 + LMT3)
    ), varList(0, 300)).solve

  println("Result:")
  println(zmap)
  val zm = zmap match {
    case Some(x) => x
    case None => throw new RuntimeException("NO MATCH")
  }
  println(zm(EBLS[BigInt, EqCoeff](z(201))))
  println(zm(EBLS[BigInt, EqCoeff](z(202))))

}