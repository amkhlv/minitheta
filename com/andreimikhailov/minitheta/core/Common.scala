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

trait IsMonomial

object Sum {
  def apply[S](xs:List[ER[S]]) :ER[S] = xs match {    
    case List()  => ZeroOfRing[S]()
    case x::rest => x + apply[S](rest)
  }
}

object Prod {
  def apply[S](xs:List[ER[S]]) :ER[S] = xs match {
    case List()  => UnitOfRing[S](1)
    case x::rest => x * apply[S](rest)
  }
}

class Rational(n: BigInt, d: BigInt)  
// This is modified from the book of Odersky, Spoon and Venners
{
  require(d != 0)
  private val g = gcd(n.abs, d.abs) * sign(d)
  val numer = n / g
  val denom = d / g
  assert(denom > 0)
  def this(n: BigInt) = this(n, 1)
  def + (that: Rational): Rational =
    new Rational(
      numer * that.denom + that.numer * denom,
      denom * that.denom
    )
  def + (i: BigInt): Rational =
    new Rational(numer + i * denom, denom)
  def - (that: Rational): Rational =
    new Rational(
      numer * that.denom - that.numer * denom,
      denom * that.denom
    )
  def - (i: BigInt): Rational =
    new Rational(numer - i * denom, denom)
  def * (that: Rational): Rational =
    new Rational(numer * that.numer, denom * that.denom)
  def * (i: BigInt): Rational =
    new Rational(numer * i, denom)
  def / (that: Rational): Rational =
    new Rational(numer * that.denom, denom * that.numer)
  def / (i: BigInt): Rational =
    new Rational(numer, denom * i)
  def isInteger: Boolean = (this.denom == 1)
  def isZero: Boolean = (this.numer == 0)
  override def toString = numer +"/"+ denom
  private def sign(a:BigInt) = if (a>0) 1 else if (a == 0) 0 else -1 
  private def gcd(a: BigInt, b: BigInt): BigInt =
    if (b == 0) a else gcd(b, a % b)
}


