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

trait ELS[S,T] { // Element of Linear Space with basis T over ER[S]   
  def *: (s:ER[S]) : ELS[S,T] = 
    if (isZero(s)) ZeroV[S,T]()
    else if (isOne(s)) this
    else (s, this) match {
      case (_, ZeroV())        => ZeroV[S,T]()
      case (s1, STimesV(y,u))  => (s1*y) *: u
      case (s1, a @ EBLS(_))   => STimesV(s1,a)
      case (s1, VSum(vx @ _*)) => ((ZeroV[S,T](): ELS[S,T]) /: vx) ( _ + s1 *: _ )
  }
  private def stv(s: ER[S], v:EBLS[S,T]) = 
    if (isZero(s)) ZeroV[S,T]()
    else if (isOne(s)) v
    else STimesV[S,T](s,v)
  private def isOne  (s:ER[S]) = (s == UnitOfRing(1))
  private def isZero (s:ER[S]) = s match {
    case ZeroOfRing()  => true
    case EmbedQ(r)     => (r.numer == 0)
    case UnitOfRing(n) => (n == 0)
    case _ => throw new RuntimeException("unmatched case in isZero()")
  }
  private def sumTerms(ts: List[ELS[S,T]]) : ELS[S,T] = ts match {
    case List()        => throw new RuntimeException("sumTerms: incorrect use")
    case List(ZeroV()) => ZeroV()
    case ZeroV()::vs   => sumTerms(vs)
    case List(v)       => v 
    case _             => VSum(ts: _*)
  }
  def +  (v:ELS[S,T]): ELS[S,T]  = (this,v) match {
    case (x, ZeroV()) => x
    case (ZeroV(), x) => x
    case (VSum(bs @ _*), VSum(a, as @ _*) ) => 
      (a /: (as.toList ::: bs.toList)) (_ + _)
    case (VSum(b, bs @ _*), VSum(as @ _*) ) => 
      (b /: (as.toList ::: bs.toList)) (_ + _)
    case (x@ EBLS(_), y @ EBLS(_))    => 
      if (x == y)  stv(UnitOfRing[S](2),x) else VSum(x,y)
    case (x@ EBLS(_), STimesV(s, y))  => 
      if (x == y)  stv((UnitOfRing[S](1) + s), x) else VSum(x, STimesV(s,y))
    case (STimesV(s,y), x @ EBLS(_))  => 
      if (x == y)  stv((UnitOfRing[S](1) + s), x) else VSum(STimesV(s,y), x)
    case (STimesV(s,y), STimesV(t,z)) => 
      if (y == z)  ((s+t) *: y) else VSum(STimesV(s,y), STimesV(t,z))
    case (VSum(as @ _*), x@ EBLS(_))   => 
      sumTerms( stv(UnitOfRing[S](1) + x.multiplicityIn(as.toList), 
                    x)::(x.termsNotSimilarIn(as.toList)) )
    case (x@ EBLS(_), VSum(as @ _*))   => 
      sumTerms( stv(UnitOfRing[S](1) + x.multiplicityIn(as.toList), 
                    x)::(x.termsNotSimilarIn(as.toList)) )
    case (VSum(as @ _*), STimesV(s,u)) => 
      sumTerms( stv((s+(u.multiplicityIn(as.toList))), 
                    u)::(u.termsNotSimilarIn(as.toList)) )
    case (STimesV(s,u), VSum(as @ _*)) => 
      sumTerms( stv((s+(u.multiplicityIn(as.toList))), 
                    u)::(u.termsNotSimilarIn(as.toList)) )
    case _  => 
      { println(" *** Did not match: " + v.toString + " + " + this.toString)
        VSum(v,this)} 
  }
  def monomials : List[ELS[S,T] with IsMonomial] = this match {   
    case ZeroV()         => List()
    case x@ STimesV(_,_) => List(x)
    case x@ EBLS(_)      => List(x)
    case VSum(xs @ _*)    => xs.toList map {
        case m:IsMonomial => m 
        case _ => throw new RuntimeException("non-monomial inside a Sum")
    }
  }
} 

case class EBLS[S,T](generator:T) extends ELS[S,T] with IsMonomial
// Element of Basis of Linear Space
{
  override def toString = generator.toString
  def multiplicityIn(as: List[ELS[S,T]]) : ER[S] = {
    val similars= as.toList filter {
      case STimesV(a, EBLS(`generator`)) => true
      case EBLS(`generator`) => true
      case _ => false }
    ((ZeroOfRing[S]():ER[S]) /: similars) (
        (acc, v) => v match {
          case STimesV(s,a) => acc + s
          case EBLS(_) => acc + UnitOfRing[S](1)
          })  }
  def termsNotSimilarIn(as: List[ELS[S,T]]) : List[ELS[S,T]] = 
    as.toList filter {
      case STimesV(a, EBLS(`generator`)) => false
      case EBLS(`generator`) => false
      case _ => true
  } 
}

case class ZeroV[S,T]() extends ELS[S,T] with IsMonomial {
  override def toString = "0V"
}

case class STimesV[S,T](x: ER[S], v:EBLS[S,T]) extends ELS[S,T] 
with IsMonomial {
  require (x match {
    case ZeroOfRing()  => false
    case UnitOfRing(_) => true
    case EmbedQ(r)     => (r.numer != 0)
    case _ => false
  })
  override def toString = "(" + x.toString + ")*:" + v
}

case class VSum[S,T](vs : ELS[S,T] *)  extends ELS[S,T] {
  require (vs.toList.length > 1)
  require ( (true /: (vs.toList map {
          case ZeroV()      => false
          case _:IsMonomial => true 
          case _            => false})) (_ && _))
  override def toString = ("(" /: vs.toList) (_ + " + " + _) + " )"
}
