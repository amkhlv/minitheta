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

/** This takes lots of time
  * to run
  */
object L3T4Eqs extends App {
  println("=== " + new java.util.Date + ": calculating pre ===")
  val pre =
    unknown(1) *: word("m^l_t^t_t^t_t") +
      unknown(2) *: word("m^t_l^t_t^t_t") +
      unknown(3) *: word("m^t_t^l_t^t_t") +
      unknown(4) *: word("m^t_t^t_l^t_t") +
      unknown(5) *: word("m^t_t^t_t^l_t") +
      unknown(6) *: word("m^t_t^t_t^t_l") +
      unknown(7) *: word("l^m_t^t_t^t_t") +
      unknown(8) *: word("t^m_l^t_t^t_t") +
      unknown(9) *: word("t^m_t^l_t^t_t") +
      unknown(10) *: word("t^m_t^t_l^t_t") +
      unknown(11) *: word("t^m_t^t_t^l_t") +
      unknown(12) *: word("t^m_t^t_t^t_l") +
      unknown(13) *: word("l^t_m^t_t^t_t") +
      unknown(14) *: word("t^l_m^t_t^t_t") +
      unknown(15) *: word("t^t_m^l_t^t_t") +
      unknown(16) *: word("t^t_m^t_l^t_t") +
      unknown(17) *: word("t^t_m^t_t^l_t") +
      unknown(18) *: word("t^t_m^t_t^t_l") +
      unknown(19) *: word("l^t_t^m_t^t_t") +
      unknown(20) *: word("t^l_t^m_t^t_t") +
      unknown(21) *: word("t^t_l^m_t^t_t") +
      unknown(22) *: word("t^t_t^m_l^t_t") +
      unknown(23) *: word("t^t_t^m_t^l_t") +
      unknown(24) *: word("t^t_t^m_t^t_l") +
      unknown(25) *: word("l^t_t^t_m^t_t") +
      unknown(26) *: word("t^l_t^t_m^t_t") +
      unknown(27) *: word("t^t_l^t_m^t_t") +
      unknown(28) *: word("t^t_t^l_m^t_t") +
      unknown(29) *: word("t^t_t^t_m^l_t") +
      unknown(30) *: word("t^t_t^t_m^t_l") +
      unknown(31) *: word("l^t_t^t_t^m_t") +
      unknown(32) *: word("t^l_t^t_t^m_t") +
      unknown(33) *: word("t^t_l^t_t^m_t") +
      unknown(34) *: word("t^t_t^l_t^m_t") +
      unknown(35) *: word("t^t_t^t_l^m_t") +
      unknown(36) *: word("t^t_t^t_t^m_l") +
      unknown(37) *: word("l^t_t^t_t^t_m") +
      unknown(38) *: word("t^l_t^t_t^t_m") +
      unknown(39) *: word("t^t_l^t_t^t_m") +
      unknown(40) *: word("t^t_t^l_t^t_m") +
      unknown(41) *: word("t^t_t^t_l^t_m") +
      unknown(42) *: word("t^t_t^t_t^l_m") +
      (unknown(101) * word("l^t_").trace) *: word("m^t_t^t_t") +
      (unknown(102) * word("l^t_").trace) *: word("t^m_t^t_t") +
      (unknown(103) * word("l^t_").trace) *: word("t^t_m^t_t") +
      (unknown(104) * word("l^t_").trace) *: word("t^t_t^m_t") +
      (unknown(105) * word("l^t_").trace) *: word("t^t_t^t_m") +
      (unknown(106) * word("m^t_").trace) *: word("l^t_t^t_t") +
      (unknown(107) * word("m^t_").trace) *: word("t^l_t^t_t") +
      (unknown(108) * word("m^t_").trace) *: word("t^t_l^t_t") +
      (unknown(109) * word("m^t_").trace) *: word("t^t_t^l_t") +
      (unknown(110) * word("m^t_").trace) *: word("t^t_t^t_l") +
      (unknown(111) * word("l^t_").trace * word("m^t_").trace) *: word("t^t_t") +
      (unknown(112) * word("l^m_t^t_").trace) *: word("t^t_t") +
      (unknown(113) * word("l^t_m^t_").trace) *: word("t^t_t") +
      (unknown(114) * word("l^t_t^m_").trace) *: word("t^t_t") +
      (unknown(115) * word("l^t_t^t_").trace) *: word("m^t_t") +
      (unknown(116) * word("l^t_t^t_").trace) *: word("t^m_t") +
      (unknown(117) * word("l^t_t^t_").trace) *: word("t^t_m") +
      (unknown(118) * word("m^t_t^t_").trace) *: word("l^t_t") +
      (unknown(119) * word("m^t_t^t_").trace) *: word("t^l_t") +
      (unknown(120) * word("m^t_t^t_").trace) *: word("t^t_l") +
      (unknown(121) * word("t^t_t^t_").trace) *: word("l^m_t") +
      (unknown(122) * word("t^t_t^t_").trace) *: word("l^t_m") +
      (unknown(123) * word("t^t_t^t_").trace) *: word("m^l_t") +
      (unknown(124) * word("t^t_t^t_").trace) *: word("t^l_m") +
      (unknown(125) * word("t^t_t^t_").trace) *: word("m^t_l") +
      (unknown(126) * word("t^t_t^t_").trace) *: word("t^m_l") +
      (unknown(127) * word("l^m_t^t_t^t_").trace) *: word("t") +
      (unknown(128) * word("l^t_m^t_t^t_").trace) *: word("t") +
      (unknown(129) * word("l^t_t^m_t^t_").trace) *: word("t") +
      (unknown(130) * word("l^t_t^t_m^t_").trace) *: word("t") +
      (unknown(131) * word("l^t_t^t_t^m_").trace) *: word("t") +
      (unknown(132) * word("l^t_t^t_t^t_").trace) *: word("m") +
      (unknown(133) * word("m^t_t^t_t^t_").trace) *: word("l") +
      (unknown(134) * word("l^m_").trace) *: word("t^t_t^t_t") +
      (unknown(135) * word("l^m_").trace * word("t^t_t^t_").trace) *: word("t") +
      (unknown(136) * word("l^t_").trace * word("m^t_t^t_").trace) *: word("t") +
      (unknown(137) * word("m^t_").trace * word("l^t_t^t_").trace) *: word("t") +
      (unknown(138) * word("l^t_").trace * word("t^t_t^t_").trace) *: word("m") +
      (unknown(139) * word("m^t_").trace * word("t^t_t^t_").trace) *: word("l")

  println("===" + new java.util.Date + ": calculating QLQRpre ===")
  val QLQRpre = QLmat(QRmat(pre))
  println("===" + new java.util.Date + ": calculating QRmll   ===")
  val QRmll = QRmat(word("m^t_l^t_t^l_t") +
    unknown(200) *: word("t^l_t^t_l^t_m"))
  println("===" + new java.util.Date + ": calculating QLlmm   ===")
  val QLlmm = QLmat(word("l^t_m^t_t^m_t") +
    unknown(200) *: word("t^m_t^t_m^t_l"))
  println("===" + new java.util.Date + ": calculating for symmetric ===")
  val zmapSymm =
    Equations(collectEquationsFromMat(
      QRmll + QLlmm + QLQRpre
    ),
      varList(0, 300)).solve
  println("===" + new java.util.Date + ": calculating for antisymmetric ===")
  val zmapAntiSymm =
    Equations(collectEquationsFromMat(
      QRmll + matScalarESCA(-1) * QLlmm + QLQRpre
    ),
      varList(0, 300)).solve

  println(new java.util.Date)
  println("symmetric solutions = " + zmapSymm)
  println("antisymmetric solutions = " + zmapAntiSymm)
  //Output: zmapSymm == None,
  // zmapAntiSymm has unknown(200) = 1
}