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

object L3T2Eqs extends App {

  /** Can we take scalar $\lambda\theta\mu\theta$, multiply it by $\lambda - \mu$
    * and get a nontrivial cohomology class?
    */
  val suspect = ((word("l^t_m^t_").trace) *: word("l")) +
    ((scalarESCA(-1) * (word("l^t_m^t_").trace)) *: word("m"))

  /** Can we get a nontrivial cohomology class from $\Phi_{up/down/mixed}$ ?
    */
  val suspect1 = word("m^t_l^t_m") + (ratioESCA(1, 4) * word("m^m_").trace) *: word("t^l_t") +
    ((ratioESCA(-1, 4) * (word("m^m_").trace))) *: (word("t^t_m") + word("m^t_t"))

  println("=== verifying that these expressions are BRST closed ===")
  assert((QLmat(suspect) + QRmat(suspect)) == mscalar(0, 1, 4))
  assert((QLmat(suspect1) + QRmat(suspect1)) == mscalar(0, 1, 4))
  println("=== OK, continuing ===")

  val qinvLL =
    unknown(101) *: word("l^l_t^t_t") +
      unknown(102) *: word("l^t_l^t_t") +
      unknown(103) *: word("l^t_t^l_t") +
      unknown(104) *: word("l^t_t^t_l") +
      unknown(105) *: word("t^l_l^t_t") +
      unknown(106) *: word("t^l_t^l_t") +
      unknown(107) *: word("t^l_t^t_l") +
      unknown(108) *: word("t^t_l^l_t") +
      unknown(109) *: word("t^t_l^t_l") +
      unknown(110) *: word("t^t_t^l_l") +
      (unknown(111) * (word("l^l_").trace)) *: word("t^t_t") +
      (unknown(112) * (word("l^t_").trace)) *: word("l^t_t") +
      (unknown(113) * (word("l^t_").trace)) *: word("t^l_t") +
      (unknown(114) * (word("l^t_").trace)) *: word("t^t_l") +
      (unknown(115) * (word("l^l_t^t_").trace)) *: word("t") +
      (unknown(116) * (word("l^t_l^t_").trace)) *: word("t") +
      (unknown(117) * (word("l^t_t^l_").trace)) *: word("t") +
      (unknown(118) * (word("l^t_t^t_").trace)) *: word("l") +
      (unknown(119) * (word("l^t_").trace) * (word("l^t_").trace)) *: word("t")

  val qinvMM =
    unknown(151) *: word("m^m_t^t_t") +
      unknown(152) *: word("m^t_m^t_t") +
      unknown(153) *: word("m^t_t^m_t") +
      unknown(154) *: word("m^t_t^t_m") +
      unknown(155) *: word("t^m_m^t_t") +
      unknown(156) *: word("t^m_t^m_t") +
      unknown(157) *: word("t^m_t^t_m") +
      unknown(158) *: word("t^t_m^m_t") +
      unknown(159) *: word("t^t_m^t_m") +
      unknown(160) *: word("t^t_t^m_m") +
      (unknown(161) * (word("m^m_").trace)) *: word("t^t_t") +
      (unknown(162) * (word("m^t_").trace)) *: word("m^t_t") +
      (unknown(163) * (word("m^t_").trace)) *: word("t^m_t") +
      (unknown(164) * (word("m^t_").trace)) *: word("t^t_m") +
      (unknown(165) * (word("m^m_t^t_").trace)) *: word("t") +
      (unknown(166) * (word("m^t_m^t_").trace)) *: word("t") +
      (unknown(167) * (word("m^t_t^m_").trace)) *: word("t") +
      (unknown(168) * (word("m^t_t^t_").trace)) *: word("m") +
      (unknown(169) * (word("m^t_").trace) * (word("m^t_").trace)) *: word("t")

  val qinvLM =
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
      (unknown(11) * (word("l^m_").trace)) *: word("t^t_t") +
      (unknown(12) * (word("l^t_").trace)) *: word("m^t_t") +
      (unknown(13) * (word("l^t_").trace)) *: word("t^m_t") +
      (unknown(14) * (word("l^t_").trace)) *: word("t^t_m") +
      (unknown(15) * (word("l^m_t^t_").trace)) *: word("t") +
      (unknown(16) * (word("l^t_m^t_").trace)) *: word("t") +
      (unknown(17) * (word("l^t_t^m_").trace)) *: word("t") +
      (unknown(18) * (word("l^t_t^t_").trace)) *: word("m") +
      (unknown(19) * (word("l^t_").trace) * (word("m^t_").trace)) *: word("t")

  val qinvML =
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
      (unknown(61) * (word("m^l_").trace)) *: word("t^t_t") +
      (unknown(62) * (word("m^t_").trace)) *: word("l^t_t") +
      (unknown(63) * (word("m^t_").trace)) *: word("t^l_t") +
      (unknown(64) * (word("m^t_").trace)) *: word("t^t_l") +
      (unknown(65) * (word("m^l_t^t_").trace)) *: word("t") +
      (unknown(66) * (word("m^t_l^t_").trace)) *: word("t") +
      (unknown(67) * (word("m^t_t^l_").trace)) *: word("t") +
      (unknown(68) * (word("m^t_t^t_").trace)) *: word("l") +
      (unknown(69) * (word("m^t_").trace) * (word("l^t_").trace)) *: word("t")

  val qinv = qinvLM + qinvML + qinvLL + qinvMM
  val zmap =
    Equations(collectEquationsFromMat(
      suspect + QLmat(qinv) + QRmat(qinv)
    ),
      varList(0, 300)).solve
  val zmap1 =
    Equations(collectEquationsFromMat(
      suspect1 + QLmat(qinv) + QRmat(qinv)
    ),
      varList(0, 300)).solve
  println(zmap)
  println("The ``Some expression'' above shows that ``suspect'' is BRST trivial")
  println(zmap1)
  println("The ``None'' above shows that ``suspect1'' is not BRST trivial")

}
