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

case class Mat[S](elemss: List[List[ER[S]]]) {
  val nrows = elemss.length
  val ncols = elemss.head.length
  //all rows should have the same length:
  require(
    (true /: elemss.tail)
      ((x: Boolean, y: List[ER[S]]) => (x && (ncols == y.length)))
  )

  def transpose: Mat[S] = Mat[S](
    for (i <- List.range(0, ncols))
    yield
      for (j <- List.range(0, nrows))
      yield elemss(j)(i)
  )

  def *(that: Mat[S]) = {
    require(ncols == that.nrows)
    Mat[S](
      for (i <- List.range(0, nrows))
      yield
        for (j <- List.range(0, that.ncols))
        yield
          Sum(for (pair <- elemss(i) zip that.transpose.elemss(j))
          yield (pair._1 * pair._2))
    )
  }

  def *:(that: ER[S]) = Mat[S](
    for (i <- List.range(0, nrows))
    yield
      for (j <- List.range(0, ncols))
      yield that * elemss(i)(j)
  )

  def +(that: Mat[S]) = {
    require((ncols == that.ncols) && (nrows == that.nrows))
    Mat[S](
      for (i <- List.range(0, nrows))
      yield
        for (j <- List.range(0, ncols))
        yield this.elemss(i)(j) + that.elemss(i)(j)
    )
  }

  val trace: ER[S] = {
    require(ncols == nrows)
    Sum(for (i <- List.range(0, nrows)) yield elemss(i)(i))
  }

  def symmetrize(): Mat[S] = (EmbedQ[S](new Rational(1,2)) *: this) +
    ((EmbedQ[S](new Rational(1,2)) *: this).transpose)
  def antisymmetrize(): Mat[S] = (EmbedQ[S](new Rational(1, 2)) *: this) +
    ((EmbedQ[S](new Rational(-1, 2)) *: this).transpose)
}

object Matrix {
  def unit[S](s: ER[S], zero: ER[S], size: Int) = Mat[S](
    for (i <- List.range(0, size))
    yield
      for (j <- List.range(0, size))
      yield
        if (i == j) s else zero
  )

  def symplectic[S](s: ER[S], zero: ER[S], halfsize: Int) = Mat[S](
    for (i <- List.range(0, halfsize); p <- List(0, 1))
    yield
      for (j <- List.range(0, halfsize); q <- List(0, 1))
      yield if (i != j) zero
      else (p, q) match {
        case (0, 0) => zero
        case (0, 1) => UnitOfRing[S](-1) * s
        case (1, 0) => s
        case (1, 1) => zero
      }
  )
}

