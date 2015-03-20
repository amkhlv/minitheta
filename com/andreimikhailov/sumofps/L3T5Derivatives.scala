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

/**
 * Created with IntelliJ IDEA.
 * User: andrei
 * Date: 29/11/12
 * Time: 12:46
 * To change this template use File | Settings | File Templates.
 */

import com.andreimikhailov.minitheta.core._
import com.andreimikhailov.brst._


object L3T5Derivatives extends App {

  def twoDerivatives(f: ESCA[BigInt, Thetas, Bosons]) = cap * Mat[ESCA[BigInt, Thetas, Bosons]](
    for (i <- List.range(0, dim)) yield for (j <- List.range(0, dim)) yield
      (scalarESCA(0) /: (
        (
          for (k <- List.range(0, dim)) yield for (l <- List.range(0, dim)) yield
            (Th(j, k) D_: (Th(i, l) D_: f)) * (cap.elemss(k)(l))
          ) flatten
        ))(_ + _)
  ) * cap

  println(new java.util.Date)
  println("====================================================================================")
  println("Calculating l3t5")
  println("1. MMMT5")
  val MMMT5 = scalarESCA(-1) * (TTM * cap * (TTM.transpose) * cup * TMup).trace
  val partial2MMM = twoDerivatives(MMMT5)
  assert(partial2MMM.antisymmetrize() == mscalar(0, 1, dim))
  // assert(partial2MMM.transpose + scalarESCA(-1) *: partial2MMM == mscalar(0, 1, dim))
  assert(QRmat(partial2MMM) == mscalar(0, 1, dim))
  //  println(partial2MMM)

  println("2. LLLT5")
  val LLLT5 = (TTL * cap * (TTL.transpose) * cup * TLup).trace
  val partial2LLL = twoDerivatives(LLLT5)
  assert(partial2LLL.antisymmetrize() == mscalar(0, 1, dim))

  println("3. LLMT5 and LMMT5")
  val LLMT5 = ratioESCA(-5, 6) * (
    ((TTL * cap * (TTL.transpose) * cup * TMup).trace) +
      ((TTL * cap * (TTM.transpose) * cup * TLup).trace) +
      ((TTM * cap * (TTL.transpose) * cup * TLup).trace) +
      (scalarESCA(-1) * (TTM * cap * (TTM.transpose) * cup * TLup).trace) +
      (scalarESCA(-1) * (TTM * cap * (TTL.transpose) * cup * TMup).trace) +
      (scalarESCA(-1) * (TTL * cap * (TTM.transpose) * cup * TMup).trace)
    )
  val LMMT5 = ratioESCA(-5, 6) * (
    ((TTL * cap * (TMdn.transpose) * (TTL.transpose) * cup).trace) +
      ((TTM * cap * (TLdn.transpose) * (TTL.transpose) * cup).trace) +
      ((TTL * cap * (TLdn.transpose) * (TTM.transpose) * cup).trace) +
      (scalarESCA(-1) * ((TTM * cap * (TLdn.transpose) * (TTM.transpose) * cup).trace)) +
      (scalarESCA(-1) * ((TTL * cap * (TMdn.transpose) * (TTM.transpose) * cup).trace)) +
      (scalarESCA(-1) * ((TTM * cap * (TMdn.transpose) * (TTL.transpose) * cup).trace))
    )

  val l3t5 = LLLT5 + LLMT5 + LMMT5 + MMMT5

  println("Verifying assertions")
  assert(QL(LLLT5) == scalarESCA(0))
  assert(QL(l3t5) + QR(l3t5) == scalarESCA(0))

  println("Calculating second derivative of l3t5")
  val l3t5partialTheta2 = twoDerivatives(l3t5)

  println("Verifying assertions")
  assert(l3t5partialTheta2.transpose + (scalarESCA(-1) *: l3t5partialTheta2) == mscalar(0, 1, dim))
  assert(QLmat(l3t5partialTheta2) + QRmat(l3t5partialTheta2) == mscalar(0, 1, dim))
  println(new java.util.Date)

  println("====================================================================================")
  println("Calculating candidate preimage")
  println("1. LLs")
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
  //println(partial2LLL + QLmat(scalarESCA(40)*:LLs(1) + scalarESCA(40)*:LLs(3) + scalarESCA(-10)*:LLs(6) + scalarESCA(-30)*:LLs(7)))

  val preimageLL = insUnknownsMat(1, LLs)

  val zmapLLL = Equations(collectEquationsFromMat(
    partial2LLL + QLmat(preimageLL)
  ),
    varList(0, 50)).solve
  zmapLLL match {
    case None => throw new RuntimeException("Even LLL part cannot be gauged away!")
    case Some(x) => println(x)
  }

  println("2. MMs")
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
  //println(partial2MMM + QRmat(scalarESCA(40)*:MMs(1) + scalarESCA(40)*:MMs(3) + scalarESCA(-10)*:MMs(6) + scalarESCA(-30)*:MMs(7)))
  //throw new RuntimeException("STOPPING")
  val preimageMM = insUnknownsMat(21, MMs)
  val zmapMMM = Equations(collectEquationsFromMat(
    partial2MMM + QRmat(preimageMM)
  ),
    varList(0, 50)).solve
  zmapMMM match {
    case None => throw new RuntimeException("Even MMM part cannot be gauged away!")
    case Some(x) => println(x)
  }
  // Map(z(28) -> (-30R)*:#, z(27) -> (-10R)*:#, z(22) -> ( + (-1R)*:z(21) + (40R)*:# ), z(24) -> ( + z(21) + (-1R)*:z(23) + (40R)*:# ))

  println("3. LMs")
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
  val preimageLM = insUnknownsMat(41, LMs)
  val preimage = preimageLL + preimageMM + preimageLM
  assert(preimage.transpose + (scalarESCA(-1) *: preimage) == mscalar(0, 1, dim))
  println(new java.util.Date)

  println("====================================================================================")
  println("Calculating our original ansatz, which is symmetrized t^l_t^t_l^m plus lower order terms ")
  println("Notice that z(220) is the coefficient of symmetrized t^l_t^t_l^m")

  println("1. Calculating L2MT3")
  val L2MT3 =
    unknown(220) *: (
      word("t^l_t^t_l^m") +
        ((ratioESCA(-1, 4) * (word("l^t_").trace)) *: word("t^l_t^m")) +
        word("m^l_t^t_l^t") +
        ((ratioESCA(-1, 4) * (word("t^l_").trace)) *: word("m^t_l^t"))
      )
  println("Verifying assertions")
  assert(QLmat(L2MT3) == mscalar(0, 1, 4))
  assert(L2MT3.transpose + (scalarESCA(-1) *: L2MT3) == mscalar(0, 1, dim))

  println("2. Calculating LM2T3")
  val LM2T3 = (
    (unknown(301) *: word("l^m_t^m_t^t")) +
      (unknown(302) *: word("l^m_t^t_m^t")) +
      (unknown(303) *: word("l^m_t^t_t^m")) +
      (unknown(304) *: word("l^t_m^t_m^t")) +
      (unknown(305) *: word("l^t_m^t_t^m")) +
      (unknown(306) *: word("l^t_t^m_t^m")) +
      (unknown(307) *: word("m^l_m^t_t^t")) +
      (unknown(308) *: word("m^l_t^m_t^t")) +
      (unknown(309) *: word("m^l_t^t_m^t")) +
      (unknown(310) *: word("m^l_t^t_t^m")) +
      (unknown(311) *: word("t^l_m^t_m^t")) +
      (unknown(312) *: word("t^l_m^t_t^m")) +
      (unknown(313) *: word("t^l_t^m_t^m")) +
      (unknown(314) *: word("m^t_l^m_t^t")) +
      (unknown(315) *: word("m^t_l^t_m^t")) +
      (unknown(316) *: word("m^t_l^t_t^m")) +
      (unknown(317) *: word("t^m_l^m_t^t")) +
      (unknown(318) *: word("t^m_l^t_m^t")) +
      (unknown(319) *: word("t^m_l^t_t^m")) +
      (unknown(320) *: word("t^t_l^m_t^m")) +
      (unknown(321) *: word("m^t_m^l_t^t")) +
      (unknown(322) *: word("m^t_t^l_m^t")) +
      (unknown(323) *: word("m^t_t^l_t^m")) +
      (unknown(324) *: word("t^m_t^l_m^t")) +
      (unknown(325) *: word("t^m_t^l_t^m")) +
      (unknown(326) *: word("t^t_m^l_m^t")) +
      (unknown(327) *: word("t^t_m^l_t^m")) +
      (unknown(328) *: word("m^t_m^t_l^t")) +
      (unknown(329) *: word("m^t_t^m_l^t")) +
      (unknown(330) *: word("t^m_t^m_l^t")) +
      (unknown(331) *: word("m^t_t^t_l^m")) +
      (unknown(332) *: word("t^m_t^t_l^m")) +
      (unknown(333) *: word("t^t_m^t_l^m")) +
      (unknown(334) *: word("t^t_t^m_l^m")) +
      (unknown(335) *: word("m^t_m^t_t^l")) +
      (unknown(336) *: word("m^t_t^m_t^l")) +
      (unknown(337) *: word("m^t_t^t_m^l")) +
      (unknown(338) *: word("t^m_t^m_t^l")) +
      (unknown(339) *: word("t^m_t^t_m^l")) +
      (unknown(340) *: word("t^t_m^t_m^l")) +
      (unknown(401) *: (word("m^m_").trace *: word("l^t_t^t"))) +
      (unknown(402) *: (word("m^m_").trace *: word("t^l_t^t"))) +
      (unknown(403) *: (word("m^m_").trace *: word("t^t_l^t"))) +
      (unknown(404) *: (word("m^m_").trace *: word("t^t_t^l"))) +
      (unknown(405) *: (word("l^m_").trace *: word("m^t_t^t"))) +
      (unknown(406) *: (word("l^m_").trace *: word("t^m_t^t"))) +
      (unknown(407) *: (word("l^m_").trace *: word("t^t_m^t"))) +
      (unknown(408) *: (word("l^m_").trace *: word("t^t_t^m"))) +
      //
      (unknown(420) *: (word("m^t_").trace *: word("l^m_t^t"))) +
      (unknown(421) *: (word("m^t_").trace *: word("l^t_m^t"))) +
      (unknown(422) *: (word("m^t_").trace *: word("l^t_t^m"))) +
      (unknown(423) *: (word("m^t_").trace *: word("t^l_m^t"))) +
      (unknown(424) *: (word("m^t_").trace *: word("t^l_t^m"))) +
      (unknown(425) *: (word("m^t_").trace *: word("t^t_l^m"))) +
      (unknown(426) *: (word("m^t_").trace *: word("m^l_t^t"))) +
      (unknown(427) *: (word("m^t_").trace *: word("m^t_l^t"))) +
      (unknown(428) *: (word("m^t_").trace *: word("m^t_t^l"))) +
      (unknown(429) *: (word("m^t_").trace *: word("t^m_l^t"))) +
      (unknown(430) *: (word("m^t_").trace *: word("t^m_t^l"))) +
      (unknown(431) *: (word("m^t_").trace *: word("t^t_m^l"))) +
      //
      (unknown(432) *: (word("l^t_").trace *: word("m^t_m^t"))) +
      (unknown(433) *: (word("l^t_").trace *: word("m^t_t^m"))) +
      (unknown(434) *: (word("l^t_").trace *: word("t^m_t^m"))) +
      //
      (unknown(450) *: (word("l^m_t^m_").trace *: word("t^t"))) +
      (unknown(452) *: (word("l^m_t^t_").trace *: word("m^t"))) +
      (unknown(453) *: (word("l^t_m^t_").trace *: word("m^t"))) +
      (unknown(454) *: (word("l^t_t^m_").trace *: word("m^t"))) +
      (unknown(455) *: (word("l^m_t^t_").trace *: word("t^m"))) +
      (unknown(456) *: (word("l^t_m^t_").trace *: word("t^m"))) +
      (unknown(457) *: (word("l^t_t^m_").trace *: word("t^m"))) +
      (unknown(458) *: (word("m^t_t^t_").trace *: word("l^m"))) +
      (unknown(459) *: (word("m^t_t^t_").trace *: word("m^l"))) +
      (unknown(481) *: ((word("l^m_").trace *
        word("m^t_").trace) *: word("t^t"))) +
      (unknown(482) *: ((word("m^m_").trace *
        word("l^t_").trace) *: word("t^t")))
    ).symmetrize()

  println("3. Calculating M3T3")
  val M3T3 = (
    unknown(491) *: (word("m^t_m^t_m^t") + word("t^m_t^m_t^m")) +
      unknown(492) *: (word("m^t_t^m_t^m") + word("m^t_m^t_t^m")) +
      unknown(493) *: ((word("m^t_").trace) *: (word("m^t_t^t") + word("t^t_t^m"))) +
      unknown(494) *: ((word("m^t_").trace) *: (word("t^m_t^t") + word("t^t_m^t"))) +
      unknown(495) *: ((word("m^t_t^t_").trace) *: (word("m^t") + scalarESCA(-1) *: word("t^m"))) +
      unknown(496) *: ((word("m^t_m^t_").trace) *: word("t^t"))
    ).symmetrize()

  assert(
    ((L2MT3 + LM2T3 + M3T3).transpose) +
      scalarESCA(-1) *: (L2MT3 + LM2T3 + M3T3) == mscalar(0, 1, dim)
  )
  println(new java.util.Date)

  println("====================================================================================")
  println("Verifying that there is something of the form L2MT3+LM2T3+M3T3 in the kernel of Q")
  val zmapToVerify = Equations(collectEquationsFromMat(
    QLmat(L2MT3 + LM2T3 + M3T3) + QRmat(L2MT3 + LM2T3 + M3T3)
  ),
    varList(0, 600)).solve
  zmapToVerify match {
    case Some(zm) => if (zm.keys.toList.contains(EBLS[BigInt, EqCoeff](z(220)))) {
      println("Value of z(220) is:")
      println(zm(EBLS[BigInt, EqCoeff](z(220))))
    } else
      println("z(220) is not determined; this is good")
    println("The keys of variable map were:")
    println(zm.keys)
    // for (k <- zm.keys.toList) println(k.toString + " ==> " +zm(k))
    case None => throw new RuntimeException("This cannot happen because all z(i) = 0 is a solution")
  }
  println(new java.util.Date)

  println("====================================================================================")
  println("Solving equations")
  val zmap =
    Equations(collectEquationsFromMat(
      l3t5partialTheta2 +
        QLmat(preimage) + QRmat(preimage) +
        L2MT3 + LM2T3 + M3T3
    ),
      varList(0, 600)).solve

  println(new java.util.Date)
  val zm = zmap match {
    case None => {
      println("Solution not found")
      Map()
    }
    case Some(x) => {
      println(
        "coefficient of symmetrized t^l_t^t_l^m + (-1/4) * (l^t_).trace * t^l_t^m = " +
          x(EBLS[BigInt, EqCoeff](z(220)))
        //coefficient of symmetrized t^l_t^t_l^m + (-1/4) * (l^t_).trace * t^l_t^m = (288R)*:#
      )
      x
    }
  }


}
