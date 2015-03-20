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

object L3T5 extends App {

  /*
   * L3T5Eqs.scala gave this:
   * scala> for (k <- zm.keys.toList ) println(k + " ==== " + zm(k))
   *
   * z(12) ==== ( + -2R|s*v|z(28) + -4R|s*v|z(26) + -1R|s*v|z(1) + -1/2R|s*v|# )
   * z(14) ==== ( + 7/3R|s*v|# + 4R|s*v|z(26) + -1R|s*v|z(4) + 2R|s*v|z(28) + 1R|s*v|z(6) )
   * z(5)  ==== ( + -4/1R|s*v|z(26) + -4R|s*v|z(28) )
   * z(13) ==== ( + -7/3R|s*v|# + -1/1R|s*v|z(6) )
   * z(19) ==== ( + -7/6R|s*v|# + -1/1R|s*v|z(25) + 2/1R|s*v|z(28) )
   * z(17) ==== ( + 4R|s*v|z(28) + 16R|s*v|z(27) )
   * z(8)  ==== ( + -1/3R|s*v|# )
   * z(23) ==== -1/12R|s*v|#
   * z(11) ==== ( + -5/2R|s*v|# + 1R|s*v|z(1) + 1R|s*v|z(4) )
   * z(22) ==== 1R|s*v|z(24)
   * z(16) ==== 2/1R|s*v|z(28)
   * z(9)  ==== ( + 7/3R|s*v|# + -1R|s*v|z(4) + 4R|s*v|z(26) )
   * z(15) ==== ( + 1/1R|s*v|z(7) + -4/1R|s*v|z(26) + 1/1R|s*v|z(4) + 1/6R|s*v|# + -4R|s*v|z(28) )
   * z(20) ==== ( + -2/1R|s*v|z(28) + -4/1R|s*v|z(27) )
   * z(2)  ==== ( + -1R|s*v|z(6) + -7/3R|s*v|# + 1R|s*v|z(4) )
   * z(3)  ==== ( + -1R|s*v|z(7) + -2/3R|s*v|# + -1R|s*v|z(4) )
   * z(10) ==== ( + -1R|s*v|z(7) + -1R|s*v|z(4) + -7/6R|s*v|# + -1R|s*v|z(1) )
  */

  println("VERIFYING CLASS:")
  val classFound =
    word("t^m_t^t_m^t_m^t_").trace +
      ratioESCA(-1, 4) * (word("t^m_t^t_m^t_").trace * word("m^t_").trace) +
      scalarESCA(-1) * word("t^l_t^t_l^t_l^t_").trace +
      (scalarESCA(-1) * ratioESCA(-1, 4)) * (word("t^l_t^t_l^t_").trace * word("l^t_").trace) +
      ratioESCA(-7, 3) * word("l^m_t^t_m^t_t^t_").trace +
      ratioESCA(7, 3) * word("m^l_t^t_l^t_t^t_").trace +
      ratioESCA(-2, 3) * word("l^m_t^t_t^m_t^t_").trace +
      ratioESCA(2, 3) * word("m^l_t^t_t^l_t^t_").trace +
      ratioESCA(-1, 3) * word("l^t_m^t_t^t_m^t_").trace +
      ratioESCA(1, 3) * word("m^t_l^t_t^t_l^t_").trace +
      ratioESCA(7, 3) * word("l^t_m^t_t^t_t^m_").trace +
      ratioESCA(-7, 3) * word("m^t_l^t_t^t_t^l_").trace +
      ratioESCA(-7, 6) * word("l^t_t^m_t^m_t^t_").trace +
      ratioESCA(7, 6) * word("m^t_t^l_t^l_t^t_").trace +
      ratioESCA(-5, 2) * word("l^t_t^m_t^t_m^t_").trace +
      ratioESCA(5, 2) * word("m^t_t^l_t^t_l^t_").trace +
      ratioESCA(-1, 2) * word("l^t_t^m_t^t_t^m_").trace +
      ratioESCA(1, 2) * word("m^t_t^l_t^t_t^l_").trace +
      ratioESCA(-7, 3) * word("l^t_t^t_m^t_m^t_").trace +
      ratioESCA(7, 3) * word("m^t_t^t_l^t_l^t_").trace +
      ratioESCA(7, 3) * word("l^t_t^t_m^t_t^m_").trace +
      ratioESCA(-7, 3) * word("m^t_t^t_l^t_t^l_").trace +
      ratioESCA(1, 6) * word("l^t_t^t_t^m_t^m_").trace +
      ratioESCA(-1, 6) * word("m^t_t^t_t^l_t^l_").trace +
      ratioESCA(-7, 6) * word("l^m_").trace * word("t^t_t^t_t^m_").trace +
      ratioESCA(7, 6) * word("m^l_").trace * word("t^t_t^t_t^l_").trace +
      ratioESCA(-1, 12) * word("l^t_").trace * word("m^t_t^m_t^t_").trace +
      ratioESCA(1, 12) * word("m^t_").trace * word("l^t_t^l_t^t_").trace

  println(QL(classFound) + QR(classFound))

  /*
   println(
       QLTcapLcupTcapTcupLcapT  +
       mscalar(-1,dim)*cupless(word("l^l"))*word("_t^t_l^t") +
       word("t^")*capless(word("l_l"))*word("^t_l^t") +
       mscalar(-1,dim)*word("t^l_t^")*capless(word("l_l"))*word("^t") +
       word("t^l_t^t_")*cupless(word("l^l"))
       )
       */
}


