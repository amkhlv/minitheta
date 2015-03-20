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

package com.andreimikhailov.minitheta.core

/** Helper object
  * T is anti-commuting generators and U commuting generators
  */
object GrAlg {
  def prodX[U](x: Map[U, Int], y: Map[U, Int]): Map[U, Int] = {
    val keys = x.keySet ++ y.keySet
    Map(
      (for (k <- keys.toList)
      yield
        (if (x.keySet.contains(k))
          (if (y.keySet.contains(k)) (k, x(k) + y(k))
          else (k, x(k)))
        else (k, y(k)))): _*) filter {
      p => (p._2 != 0)
    }
  }

  def signTh[T <% Ordered[T]](x: List[T], y: List[T]): Int = (x, y) match {
    case (List(), _) => 1
    case (_, List()) => 1
    case (a :: rest, bs) =>
      if (bs.filter(_ == a) == List())
        ((1 /: (bs filter (_ < a)))((u, v) => u * (-1))) * signTh[T](rest, bs)
      else 0
  }
}

/** Element of FREE Super Commutative Algebra 
  * (I should have called it EFSCA?)
  */
case class ESCA[S, T <% Ordered[T], U]
(vec: ELS[S, Tuple2[Set[T], Map[U, Int]]])
  extends ER[ESCA[S, T, U]] {
  type MONOMIALS = Tuple2[Set[T], Map[U, Int]]

  implicit def ELStoESCA(x: ELS[S, MONOMIALS]): ESCA[S, T, U] = ESCA[S, T, U](x)

  implicit def ESCAtoELS(x: ESCA[S, T, U]): ELS[S, Tuple2[Set[T], Map[U, Int]]] =
    x match {
      case ESCA(y) => y
    }

  override def incr(that: BigInt) = if (that == 0) this
  else
    this + ESCA[S, T, U](
      STimesV(UnitOfRing(that), EBLS[S, MONOMIALS]((Set(), Map())))
    )


  override def incr(that: Rational) = if (that.isZero) this
  else if (that.isInteger) this.incr(that.numer)
  else
    this + ESCA[S, T, U](
      STimesV(EmbedQ(that), EBLS[S, MONOMIALS]((Set(), Map())))
    )

  override def multiply(that: BigInt) = if (that == 0) ZeroV[S, MONOMIALS]()
  else
    this * ESCA[S, T, U](STimesV(UnitOfRing(that),
      EBLS[S, MONOMIALS]((Set(), Map()))))

  override def multiply(that: Rational) = if (that.isZero) ZeroV[S, MONOMIALS]()
  else if (that.isInteger) this.multiply(that.numer)
  else
    this * ESCA[S, T, U](STimesV(EmbedQ(that),
      EBLS[S, MONOMIALS]((Set(), Map()))))

  override def *(that: ER[ESCA[S, T, U]]): ER[ESCA[S, T, U]] =
    (this, that) match {
      case (_, ZeroOfRing()) => ZeroOfRing[ESCA[S, T, U]]()
      case (z, x@ESCA(_)) => z * x
    }

  override def +(that: ER[ESCA[S, T, U]]): ER[ESCA[S, T, U]] = that match {
    case x@ESCA(_) => this.vec + x.vec
    case ZeroOfRing() => this.vec
  }

  def *(that: ESCA[S, T, U]): ESCA[S, T, U] = (this.vec, that.vec) match {
    case (ZeroV(), _) => ZeroV[S, MONOMIALS]()
    case (_, ZeroV()) => ZeroV[S, MONOMIALS]()
    case (VSum(w@_*), x) =>
      ((ZeroV[S, MONOMIALS](): ELS[S, MONOMIALS]) /: w.toList)((a, b) => a + (b * x))
    case (x, VSum(w@_*)) =>
      ((ZeroV[S, MONOMIALS](): ELS[S, MONOMIALS]) /: w.toList)((a, b) => a + (x * b))
    case (EBLS(x), EBLS(y)) => {
      val sgn = GrAlg.signTh[T](x._1.toList, y._1.toList)
      if (sgn == 0) ZeroV[S, MONOMIALS]()
      else UnitOfRing[S](sgn) *: EBLS[S, MONOMIALS]((x._1 ++ y._1,
        GrAlg.prodX(x._2, y._2)))
    }
    case (STimesV(s, v), a@EBLS(_)) => s *: (v * a)
    case (a@EBLS(_), STimesV(s, v)) => s *: (a * v)
    case (STimesV(sl, vl), STimesV(sr, vr)) => (sl * sr) *: (vl * vr)
  }

  def D_:(that: T): ESCA[S, T, U] = (that, this.vec) match {
    case (_, ZeroV()) => ZeroV[S, MONOMIALS]()
    case (x, VSum(w@_*)) =>
      ((ZeroV[S, MONOMIALS](): ELS[S, MONOMIALS]) /: w.toList)(
        (a, b) => a + (x D_: b))
    case (t, EBLS(x)) => if (x._1.contains(t)) {
      val sgn = GrAlg.signTh[T](Set(t).toList, (x._1 - t).toList)
      UnitOfRing[S](sgn) *: EBLS[S, MONOMIALS](x._1 - t, x._2)
    } else ZeroV[S, MONOMIALS]()
    case (t, STimesV(s, v)) => s *: (t D_: v)
  }
}
  
