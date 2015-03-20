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

trait ER[S] // Element of the Ring S 
{
  def + (that:ER[S]) : ER[S]  
  def * (that:ER[S]) : ER[S]
  def incr (that: BigInt): ER[S]
  def incr (that: Rational): ER[S] = 
    throw new RuntimeException("This is only a ring, not a field")
  def multiply (that: BigInt): ER[S]
  def multiply (that: Rational): ER[S] = 
    throw new RuntimeException("This is only a ring, not a field")
  def inv : ER[S] = 
    throw new RuntimeException("This is only a ring, not a field")   
} 

case class UnitOfRing[S](n:BigInt) extends ER[S] with IsMonomial {
  override def + (that:ER[S]) = that.incr(n)
  override def * (that:ER[S]) = that.multiply(n)
  override def incr (that: BigInt) = 
    if ((that + n) == 0) ZeroOfRing[S]() else UnitOfRing[S](that + n)
  override def incr (that: Rational) = {
    val p = (that + new Rational(n,1))
    if (p == new Rational(0,1)) ZeroOfRing[S]() else EmbedQ[S](p)
  }
  override def multiply (that: BigInt) = {
    val p = that * n
    if ( p == 0 ) ZeroOfRing[S]() else UnitOfRing[S](that * n)
  }
  override def multiply (that: Rational) = {
    val res = that * new Rational(n,1)
    if      (res == new Rational (0,1)) ZeroOfRing[S]()
    else if (res.denom == 1)            UnitOfRing[S](res.numer) 
    else                                EmbedQ[S](res)
  }
  override def toString = n.toString() + "R"
  override def inv = EmbedQ[S](new Rational(1,n))
}

case class EmbedQ[S](r:Rational) extends ER[S] with IsMonomial {
  //TODO: check if result is zero if so cast as ZeroOfRing[S] 
  override def + (that:ER[S]) = that.incr(r)
  override def * (that:ER[S]) = that.multiply(r)
  override def incr (that:BigInt) = {
    val res = r + new Rational(that,1) 
    if (res.numer == 0) ZeroOfRing[S]() 
    else if (res.denom == 1) UnitOfRing[S](res.numer) else EmbedQ[S](res)
  }
  override def incr (that:Rational) = {
    val res = r + that
    if (res.numer == 0) ZeroOfRing[S]() 
    else if (res.denom == 1) UnitOfRing[S](res.numer) else EmbedQ[S](res)
  }
  override def multiply(that:BigInt) = {
    val res = r * new Rational(that,1)
    if (res.numer == 0) ZeroOfRing[S]() 
    else if (res.denom == 1) UnitOfRing[S](res.numer) else EmbedQ[S](res)
  }
  override def multiply(that:Rational) = {
    val res = r * that
    if (res.numer == 0) ZeroOfRing[S]() 
    else if (res.denom == 1) UnitOfRing[S](res.numer) else EmbedQ[S](res)
  }
  override def inv = EmbedQ[S](new Rational(r.denom, r.numer))
  override def toString = r.toString() + "R"
}

case class ZeroOfRing[S]() extends ER[S] with IsMonomial {
  override def + (that:ER[S]) = that
  override def * (that:ER[S]) = ZeroOfRing[S]()
  override def incr (that: BigInt) = UnitOfRing[S](that)
  override def incr (that: Rational) = EmbedQ[S](that)
  override def multiply (that: BigInt) = ZeroOfRing[S]()
  override def multiply (that: Rational) = ZeroOfRing[S]()
  override def inv = throw new RuntimeException("Division by zero")
  override def toString = "0ofR"
}
