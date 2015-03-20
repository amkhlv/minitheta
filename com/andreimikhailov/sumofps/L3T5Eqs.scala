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

object L3T5Eqs extends App {

  println("=== " + new java.util.Date + ": calculating LM2T5 ===")
  val LM2T5 =
    unknown(1) * word("l^m_t^m_t^t_t^t_").trace +
      unknown(2) * word("l^m_t^t_m^t_t^t_").trace +
      unknown(3) * word("l^m_t^t_t^m_t^t_").trace +
      unknown(4) * word("l^m_t^t_t^t_m^t_").trace +
      unknown(5) * word("l^m_t^t_t^t_t^m_").trace +
      unknown(6) * word("l^t_m^t_m^t_t^t_").trace +
      unknown(7) * word("l^t_m^t_t^m_t^t_").trace +
      unknown(8) * word("l^t_m^t_t^t_m^t_").trace +
      unknown(9) * word("l^t_m^t_t^t_t^m_").trace +
      unknown(10) * word("l^t_t^m_t^m_t^t_").trace +
      unknown(11) * word("l^t_t^m_t^t_m^t_").trace +
      unknown(12) * word("l^t_t^m_t^t_t^m_").trace +
      unknown(13) * word("l^t_t^t_m^t_m^t_").trace +
      unknown(14) * word("l^t_t^t_m^t_t^m_").trace +
      unknown(15) * word("l^t_t^t_t^m_t^m_").trace +
      unknown(16) * word("l^t_t^t_t^m_").trace * word("t^m_").trace +
      unknown(17) * word("l^t_t^m_").trace * word("t^t_t^m_").trace +
      unknown(18) * word("l^t_m^t_").trace * word("t^t_t^m_").trace +
      unknown(19) * word("l^m_").trace * word("t^t_t^t_t^m_").trace +
      unknown(20) * word("l^m_t^m_").trace * word("t^t_t^t_").trace +
      unknown(21) * word("l^t_t^t_").trace * word("m^t_m^t_").trace +
      unknown(22) * word("l^t_").trace * word("m^t_m^t_t^t_").trace +
      unknown(23) * word("l^t_").trace * word("m^t_t^m_t^t_").trace +
      unknown(24) * word("l^t_").trace * word("m^t_t^t_m^t_").trace +
      unknown(25) * word("l^m_").trace * word("m^t_t^t_t^t_").trace +
      unknown(26) * word("l^t_t^t_t^t_").trace * word("m^m_").trace +
      unknown(27) * word("l^t_").trace * word("m^m_").trace * word("t^t_t^t_").trace +
      unknown(28) * word("l^m_").trace * word("m^t_").trace * word("t^t_t^t_").trace

  println("=== " + new java.util.Date + ": calculating L2MT5 ===")
  val L2MT5 =
    unknown(1) * word("m^l_t^l_t^t_t^t_").trace +
      unknown(2) * word("m^l_t^t_l^t_t^t_").trace +
      unknown(3) * word("m^l_t^t_t^l_t^t_").trace +
      unknown(4) * word("m^l_t^t_t^t_l^t_").trace +
      unknown(5) * word("m^l_t^t_t^t_t^l_").trace +
      unknown(6) * word("m^t_l^t_l^t_t^t_").trace +
      unknown(7) * word("m^t_l^t_t^l_t^t_").trace +
      unknown(8) * word("m^t_l^t_t^t_l^t_").trace +
      unknown(9) * word("m^t_l^t_t^t_t^l_").trace +
      unknown(10) * word("m^t_t^l_t^l_t^t_").trace +
      unknown(11) * word("m^t_t^l_t^t_l^t_").trace +
      unknown(12) * word("m^t_t^l_t^t_t^l_").trace +
      unknown(13) * word("m^t_t^t_l^t_l^t_").trace +
      unknown(14) * word("m^t_t^t_l^t_t^l_").trace +
      unknown(15) * word("m^t_t^t_t^l_t^l_").trace +
      unknown(16) * word("m^t_t^t_t^l_").trace * word("t^l_").trace +
      unknown(17) * word("m^t_t^l_").trace * word("t^t_t^l_").trace +
      unknown(18) * word("m^t_l^t_").trace * word("t^t_t^l_").trace +
      unknown(19) * word("m^l_").trace * word("t^t_t^t_t^l_").trace +
      unknown(20) * word("m^l_t^l_").trace * word("t^t_t^t_").trace +
      unknown(21) * word("m^t_t^t_").trace * word("l^t_l^t_").trace +
      unknown(22) * word("m^t_").trace * word("l^t_l^t_t^t_").trace +
      unknown(23) * word("m^t_").trace * word("l^t_t^l_t^t_").trace +
      unknown(24) * word("m^t_").trace * word("l^t_t^t_l^t_").trace +
      unknown(25) * word("m^l_").trace * word("l^t_t^t_t^t_").trace +
      unknown(26) * word("m^t_t^t_t^t_").trace * word("l^l_").trace +
      unknown(27) * word("m^t_").trace * word("l^l_").trace * word("t^t_t^t_").trace +
      unknown(28) * word("m^l_").trace * word("l^t_").trace * word("t^t_t^t_").trace

  println("=== " + new java.util.Date + ": solving for coefficients ===")
  /*
  // This should give None:
  val zmap =  
      Equations(collectEquations(
      QL(            word("t^m_t^t_m^t_m^t_").trace +           ratioESCA(-1,4)*(word("t^m_t^t_m^t_").trace * word("m^t_").trace)) +
      QR( LM2T5 ) + QL( LM2T5 ) + QR( L2MT5 ) + QL( L2MT5 ) +
      QR(            word("t^l_t^t_l^t_l^t_").trace +           ratioESCA(-1,4)*(word("t^l_t^t_l^t_").trace * word("l^t_").trace))
      ), for (i <- List.range(0,300)) yield EBLS[BigInt,EqCoeff](z(i))).solve
    */

  // This should give Some:
  val zmap =
    Equations(collectEquations(
      QL(word("t^m_t^t_m^t_m^t_").trace + ratioESCA(-1, 4) * (word("t^m_t^t_m^t_").trace * word("m^t_").trace)) +
        QR(LM2T5) + QL(LM2T5) + QR(scalarESCA(-1) * L2MT5) + QL(scalarESCA(-1) * L2MT5) +
        QR(scalarESCA(-1) * word("t^l_t^t_l^t_l^t_").trace + (scalarESCA(-1) * ratioESCA(-1, 4)) * (word("t^l_t^t_l^t_").trace * word("l^t_").trace))
    ), varList(0, 300)).solve

  println(new java.util.Date)
  println("zmap = " + zmap.toString)
  val zm = zmap match {
    case Some(x) => x
    case None => throw new RuntimeException("NO MATCH")
  }
  for (k <- zm.keys) println(k + "--->  " + zm(k))

}