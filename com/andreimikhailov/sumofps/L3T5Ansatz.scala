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
 * Date: 21/11/12
 * Time: 14:52
 * To change this template use File | Settings | File Templates.
 */

import com.andreimikhailov.minitheta.core._
import com.andreimikhailov.brst._

object L3T5Ansatz extends App {

  /** This is the ansatz for $\lambda^3\theta^5$
    */
  val ansatz =
    ((TTL * cap * (TTL.transpose) * cup * TLup).trace) +
      unknown(1) * (
        ((TTL * cap * (TTL.transpose) * cup * TMup).trace) +
          ((TTL * cap * (TTM.transpose) * cup * TLup).trace) +
          ((TTM * cap * (TTL.transpose) * cup * TLup).trace) +
          (scalarESCA(-1) * (TTM * cap * (TTM.transpose) * cup * TLup).trace) +
          (scalarESCA(-1) * (TTM * cap * (TTL.transpose) * cup * TMup).trace) +
          (scalarESCA(-1) * (TTL * cap * (TTM.transpose) * cup * TMup).trace)
        ) +
      unknown(2) * (
        ((TTL * cap * (TMdn.transpose) * (TTL.transpose) * cup).trace) +
          ((TTM * cap * (TLdn.transpose) * (TTL.transpose) * cup).trace) +
          ((TTL * cap * (TLdn.transpose) * (TTM.transpose) * cup).trace) +
          (scalarESCA(-1) * ((TTM * cap * (TLdn.transpose) * (TTM.transpose) * cup).trace)) +
          (scalarESCA(-1) * ((TTL * cap * (TMdn.transpose) * (TTM.transpose) * cup).trace)) +
          (scalarESCA(-1) * ((TTM * cap * (TMdn.transpose) * (TTL.transpose) * cup).trace))
        ) +
      (scalarESCA(-1) * (TTM * cap * (TTM.transpose) * cup * TMup).trace)
  val zmap =
    Equations(collectEquations(
      QL(ansatz) + QR(ansatz)
    ), varList(0, 100)).solve
  println(zmap)
  // This gives unknown(1) and unknown(2) both equal -5/6
}
