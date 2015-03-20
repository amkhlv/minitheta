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


object L3T3Eqs extends App {

  println("=== verifying assertion ===")
  assert(QLmat(traceless(word("l^t_")) * word("t^l_t")) == mscalar(0, 1, 4))

  val L2MT3 = cupless(
    word("t^l_t^t_l^m") +
      ((ratioESCA(-1, 4) * (word("l^t_").trace)) *: word("t^l_t^m")) +
      word("m^l_t^t_l^t") +
      ((ratioESCA(-1, 4) * (word("t^l_").trace)) *: word("m^t_l^t"))
  )

  println("=== verifying assertion ===")
  assert(QLmat(L2MT3) == mscalar(0, 1, 4))

  println("=== calculating L3T3 ===")

  val LM2T3 = cupless(
    (unknown(1) *: word("l^m_t^m_t^t")) +
      (unknown(2) *: word("l^m_t^t_m^t")) +
      (unknown(3) *: word("l^m_t^t_t^m")) +
      (unknown(4) *: word("l^t_m^t_m^t")) +
      (unknown(5) *: word("l^t_m^t_t^m")) +
      (unknown(6) *: word("l^t_t^m_t^m")) +
      (unknown(7) *: word("m^l_m^t_t^t")) +
      (unknown(8) *: word("m^l_t^m_t^t")) +
      (unknown(9) *: word("m^l_t^t_m^t")) +
      (unknown(10) *: word("m^l_t^t_t^m")) +
      (unknown(11) *: word("t^l_m^t_m^t")) +
      (unknown(12) *: word("t^l_m^t_t^m")) +
      (unknown(13) *: word("t^l_t^m_t^m")) +
      (unknown(14) *: word("m^t_l^m_t^t")) +
      (unknown(15) *: word("m^t_l^t_m^t")) +
      (unknown(16) *: word("m^t_l^t_t^m")) +
      (unknown(17) *: word("t^m_l^m_t^t")) +
      (unknown(18) *: word("t^m_l^t_m^t")) +
      (unknown(19) *: word("t^m_l^t_t^m")) +
      (unknown(20) *: word("t^t_l^m_t^m")) +
      (unknown(21) *: word("m^t_m^l_t^t")) +
      (unknown(22) *: word("m^t_t^l_m^t")) +
      (unknown(23) *: word("m^t_t^l_t^m")) +
      (unknown(24) *: word("t^m_t^l_m^t")) +
      (unknown(25) *: word("t^m_t^l_t^m")) +
      (unknown(26) *: word("t^t_m^l_m^t")) +
      (unknown(27) *: word("t^t_m^l_t^m")) +
      (unknown(28) *: word("m^t_m^t_l^t")) +
      (unknown(29) *: word("m^t_t^m_l^t")) +
      (unknown(30) *: word("t^m_t^m_l^t")) +
      (unknown(31) *: word("m^t_t^t_l^m")) +
      (unknown(32) *: word("t^m_t^t_l^m")) +
      (unknown(33) *: word("t^t_m^t_l^m")) +
      (unknown(34) *: word("t^t_t^m_l^m")) +
      (unknown(35) *: word("m^t_m^t_t^l")) +
      (unknown(36) *: word("m^t_t^m_t^l")) +
      (unknown(37) *: word("m^t_t^t_m^l")) +
      (unknown(38) *: word("t^m_t^m_t^l")) +
      (unknown(39) *: word("t^m_t^t_m^l")) +
      (unknown(40) *: word("t^t_m^t_m^l")) +
      (unknown(101) *: (word("m^m_").trace *: word("l^t_t^t"))) +
      (unknown(102) *: (word("m^m_").trace *: word("t^l_t^t"))) +
      (unknown(103) *: (word("m^m_").trace *: word("t^t_l^t"))) +
      (unknown(104) *: (word("m^m_").trace *: word("t^t_t^l"))) +
      (unknown(105) *: (word("l^m_").trace *: word("m^t_t^t"))) +
      (unknown(106) *: (word("l^m_").trace *: word("t^m_t^t"))) +
      (unknown(107) *: (word("l^m_").trace *: word("t^t_m^t"))) +
      (unknown(108) *: (word("l^m_").trace *: word("t^t_t^m"))) +
      //
      (unknown(120) *: (word("m^t_").trace *: word("l^m_t^t"))) +
      (unknown(121) *: (word("m^t_").trace *: word("l^t_m^t"))) +
      (unknown(122) *: (word("m^t_").trace *: word("l^t_t^m"))) +
      (unknown(123) *: (word("m^t_").trace *: word("t^l_m^t"))) +
      (unknown(124) *: (word("m^t_").trace *: word("t^l_t^m"))) +
      (unknown(125) *: (word("m^t_").trace *: word("t^t_l^m"))) +
      (unknown(126) *: (word("m^t_").trace *: word("m^l_t^t"))) +
      (unknown(127) *: (word("m^t_").trace *: word("m^t_l^t"))) +
      (unknown(128) *: (word("m^t_").trace *: word("m^t_t^l"))) +
      (unknown(129) *: (word("m^t_").trace *: word("t^m_l^t"))) +
      (unknown(130) *: (word("m^t_").trace *: word("t^m_t^l"))) +
      (unknown(131) *: (word("m^t_").trace *: word("t^t_m^l"))) +
      //
      (unknown(132) *: (word("l^t_").trace *: word("m^t_m^t"))) +
      (unknown(133) *: (word("l^t_").trace *: word("m^t_t^m"))) +
      (unknown(134) *: (word("l^t_").trace *: word("t^m_t^m"))) +
      //
      (unknown(150) *: (word("l^m_t^m_").trace *: word("t^t"))) +
      (unknown(152) *: (word("l^m_t^t_").trace *: word("m^t"))) +
      (unknown(153) *: (word("l^t_m^t_").trace *: word("m^t"))) +
      (unknown(154) *: (word("l^t_t^m_").trace *: word("m^t"))) +
      (unknown(155) *: (word("l^m_t^t_").trace *: word("t^m"))) +
      (unknown(156) *: (word("l^t_m^t_").trace *: word("t^m"))) +
      (unknown(157) *: (word("l^t_t^m_").trace *: word("t^m"))) +
      (unknown(200) *: ((word("l^m_").trace *
        word("m^t_").trace) *: word("t^t"))) +
      (unknown(201) *: ((word("m^m_").trace *
        word("l^t_").trace) *: word("t^t")))
  )

  val M3T3 = cupless(
    unknown(251) *: (word("m^t_m^t_m^t") + word("t^m_t^m_t^m")) +
      unknown(252) *: (word("m^t_t^m_t^m") + word("m^t_m^t_t^m")) +
      unknown(253) *: ((word("m^t_").trace) *: (word("m^t_t^t") + word("t^t_t^m"))) +
      unknown(254) *: ((word("m^t_").trace) *: (word("t^m_t^t") + word("t^t_m^t"))) +
      unknown(255) *: ((word("m^t_t^t_").trace) *: (word("m^t") + scalarESCA(-1) *: word("t^m"))) +
      unknown(256) *: ((word("m^t_m^t_").trace) *: word("t^t"))
  )

  println("Calculating candidate preimage")
  println("LLs")
  val LLs = List(
    word("l^t_l^t_t^t") + scalarESCA(-1) *: word("t^t_t^l_t^l"),
    word("l^t_t^l_t^t") + scalarESCA(-1) *: word("t^t_l^t_t^l"),
    word("l^t_t^t_l^t") + scalarESCA(-1) *: word("t^l_t^t_t^l"),
    word("t^l_t^l_t^t") + scalarESCA(-1) *: word("t^t_l^t_l^t"),
    (word("l^t_l^t_").trace) *: word("t^t"),
    (word("l^t_t^t_").trace) *: (word("l^t") + scalarESCA(-1) *: word("t^l")),
    (word("l^t_").trace) *: (word("l^t_t^t") + word("t^t_t^l")),
    (word("l^t_").trace) *: (word("t^l_t^t") + word("t^t_l^t"))
  )
  val preimageLL = insUnknownsMat(301, LLs)


  println("MMs")
  val MMs = List(
    word("m^t_m^t_t^t") + scalarESCA(-1) *: word("t^t_t^m_t^m"),
    word("m^t_t^m_t^t") + scalarESCA(-1) *: word("t^t_m^t_t^m"),
    word("m^t_t^t_m^t") + scalarESCA(-1) *: word("t^m_t^t_t^m"),
    word("t^m_t^m_t^t") + scalarESCA(-1) *: word("t^t_m^t_m^t"),
    (word("m^t_m^t_").trace) *: word("t^t"),
    (word("m^t_t^t_").trace) *: (word("m^t") + scalarESCA(-1) *: word("t^m")),
    (word("m^t_").trace) *: (word("m^t_t^t") + word("t^t_t^m")),
    (word("m^t_").trace) *: (word("t^m_t^t") + word("t^t_m^t"))
  )
  val preimageMM = insUnknownsMat(321, MMs)

  println("LMs")
  val LMs = List(
    word("m^l_t^t_t^t") + scalarESCA(-1) *: word("t^t_t^t_l^m"),
    word("m^t_l^t_t^t") + scalarESCA(-1) *: word("t^t_t^l_t^m"),
    word("m^t_t^l_t^t") + scalarESCA(-1) *: word("t^t_l^t_t^m"),
    word("m^t_t^t_l^t") + scalarESCA(-1) *: word("t^l_t^t_t^m"),
    word("m^t_t^t_t^l") + scalarESCA(-1) *: word("l^t_t^t_t^m"),
    word("t^m_l^t_t^t") + scalarESCA(-1) *: word("t^t_t^l_m^t"),
    word("t^m_t^l_t^t") + scalarESCA(-1) *: word("t^t_l^t_m^t"),
    word("t^m_t^t_l^t") + scalarESCA(-1) *: word("t^l_t^t_m^t"),
    word("t^m_t^t_t^l") + scalarESCA(-1) *: word("l^t_t^t_m^t"),
    word("t^t_m^l_t^t") + scalarESCA(-1) *: word("t^t_l^m_t^t"),
    word("t^t_m^t_l^t") + scalarESCA(-1) *: word("t^l_t^m_t^t"),
    word("t^t_m^t_t^l") + scalarESCA(-1) *: word("l^t_t^m_t^t"),
    word("t^t_t^m_l^t") + scalarESCA(-1) *: word("t^l_m^t_t^t"),
    word("t^t_t^m_t^l") + scalarESCA(-1) *: word("l^t_m^t_t^t"),
    word("t^t_t^t_m^l") + scalarESCA(-1) *: word("l^m_t^t_t^t"),
    (word("l^m_t^t_").trace) *: word("t^t"),
    (word("l^t_m^t_").trace) *: word("t^t"),
    (word("l^t_t^m_").trace) *: word("t^t"),
    (word("l^t_t^t_").trace) *: (word("m^t") + scalarESCA(-1) *: word("t^m")),
    (word("m^t_t^t_").trace) *: (word("l^t") + scalarESCA(-1) *: word("t^l")),
    ((word("l^t_").trace) * (word("m^t_").trace)) *: word("t^t"),
    (word("l^t_").trace) *: (word("m^t_t^t") + word("t^t_t^m")),
    (word("l^t_").trace) *: (word("t^m_t^t") + word("t^t_m^t")),
    (word("m^t_").trace) *: (word("l^t_t^t") + word("t^t_t^l")),
    (word("m^t_").trace) *: (word("t^l_t^t") + word("t^t_l^t")),
    (word("t^t_t^t_").trace) *: (word("l^m") + scalarESCA(-1) *: word("m^l"))
  )
  val preimageLM = insUnknownsMat(341, LMs)

  val preimage = preimageLL + preimageMM + preimageLM
  assert(preimage + scalarESCA(-1) *: preimage.transpose == mscalar(0, 1, dim))

  println("=== deriving equations ===")
  /*
  val zmap =
    Equations(collectEquationsFromMat(
      QLmat(LM2T3) + QRmat(L2MT3)
    ), for (i <- List.range(0, 300)) yield EBLS[BigInt, EqCoeff](z(i))).solve
    */
  val L3tot = L2MT3 + LM2T3 + M3T3
  val zmap =
    Equations(collectEquationsFromMat(
      L3tot + (L3tot.transpose) + QLmat(preimage) + QRmat(preimage)
    ), for (i <- List.range(0, 500)) yield EBLS[BigInt, EqCoeff](z(i))).solve

  println("")
  println("Solution:")
  println(zmap)

}