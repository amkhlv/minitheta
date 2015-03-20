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

package com.andreimikhailov.minitheta.examples

import  com.andreimikhailov.minitheta.core._

object Examples extends App {
  /** Element of basis of the linear space:
    *
    * @param s is the name of the variable ''e.g.'' "x" or "y"
    * @return  element of the basis
    */
  def V(s:String): ELS[BigInt, String] = EBLS[BigInt, String](s)

  /** Integer coefficients:
    *
    * @param m integer
    * @return  m times the unit of the ring
    */
  def N(m:BigInt): UnitOfRing[BigInt] = UnitOfRing(m)

  /** Rational coefficients:
    *
    * @param m numerator
    * @param n denominator
    * @return {{{ $m\over n$ times unit of the ring }}}
    */
  def R(m:BigInt, n:BigInt): EmbedQ[BigInt] = EmbedQ[BigInt](new Rational(m,n))

  /** Sample system of linear equation
    * Here x and y are variables, while a and b are constants.
    * The equations are:
    * {{{
    *   x + 7y + a = 0
    *   y + {1\over 13}b = 0
    * }}}
    */
  val eqs = Equations(
    List(V("x") + N(7)*:V("y") + V("a"), V("y") + R(1,13)*:V("b")),
    List("x","y") map ( u => EBLS[BigInt,String](u))
  )
  println(eqs.solve)


  /** Fermionic variables. When we define fermionic variables, we have to explain how they are ordered.
    * This is because of the anticommutativity of fermions:
    * {{{ \theta_1 \theta_2 = - \theta_2 \theta_1 }}}
    * We have to keep track of how thetas are ordered, so we have to '''define''' their order.
    * Internally, monomials correspond to sets of thetas. For example:
    * {{{ \theta_1\theta_2\theta_5 }}}
    * is a set (Th(1),Th(2),Th(5))
    */
  abstract class Thetas
  case class Th(m:Int) extends Thetas
  case class ThOrd(m:Int) extends Ordered[Thetas]{
    override def compare(that:Thetas):Int = that match {
      case Th(p) => (p - m)
    }
  }
  implicit def ordTheta(x:Thetas): Ordered[Thetas] = x match {case Th(p) => ThOrd(p)}

  /** Now introduce some bosonic (''i.e.'' usual, commutative) variables
    */
  abstract class Bosons
  case class X(m:Int) extends Bosons

  /** Super-monomial.
    * As we already explained, a monomial of fermionic variables is just a set of Thetas
    * And what is a monomial of bosonic variables?
    * It is a '''Map''', more precisely [[scala.collection.mutable.Map]];
    * Every bosonic variable is mapped to its power.
    * For example
    * {{{ $X_1^5 X_2^8 X_5$ }}}
    * corresponds to Map(X(1) -> 5, X(2) -> 8, X(5) -> 1)
    * (If the power is zero, ''i.e.'' the monomial does not contain a particular variable,
    *  then we simply drop this key.)
    */
  type BosonsFermions = Tuple2[Set[Thetas], Map[Bosons,Int]]

  /** Linear combination of super-monomials
    */
  type SUSY = ESCA[BigInt,Thetas,Bosons]

  // We need a few casting rules:
  implicit def ESCAtoELS(x: SUSY) :ELS[BigInt, BosonsFermions] = x match {case ESCA(y) => y}
  implicit def ERESCAtoESCA(x:ER[SUSY]):SUSY = x match { case y@ ESCA(_) => y }

  /** This is just a shortcut to treat an element of linear space as an element of free supercommutative algebra.
    * @param x element of the '''linear space''' of super-polynomials
    * @return  the same x, but as an element of free supercommutative algebra (a tautology)
    */
  def esca(x:ELS[BigInt,BosonsFermions]) = ESCA[BigInt,Thetas,Bosons](x)

  /** This is a shortcut to contruct a base monomial
    *
    * @param t set of fermions
    * @param m map of bosons
    * @return corresponding element of the basis of the linear space of super-polynomials
    */
  def monoBase(t:Set[Thetas], m:Map[Bosons,Int]):EBLS[BigInt, BosonsFermions] =
    EBLS[BigInt, BosonsFermions]((t,m))
  def mono(t:Set[Thetas], m:Map[Bosons,Int]):ESCA[BigInt,Thetas,Bosons] =
    esca(monoBase(t,m))

  /** This is a sample super-monomial
    */
  val M1 = mono(Set(Th(1),Th(5)), Map(X(1) -> 2))

  /** And another one
    */
  val M2 = mono(Set(Th(2),Th(6)), Map(X(1) -> 5, X(2) -> 1))

  /** And their product.
    * Notice that there is a minus sign, because we had to exchange Th(2) and Th(5)
    * to bring this to our defined ordering
    */
  val M12 = M1 * M2

  println(M1.toString + "  multiplied by  " + M2.toString + "  gives  " + M12.toString)

  /** This is the left derivative of the product M1 M2
    * with respect to Th(5).
    * Notice that the sign is now plus, because we got another minus when carrying the derivative
    * across Th(1)
    */
  val DM12 = Th(5) D_: M12

  println("Derivative of " + M12.toString + " with respect to " + Th(5).toString + " gives " + DM12.toString)


}