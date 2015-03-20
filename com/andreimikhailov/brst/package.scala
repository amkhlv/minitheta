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

package com.andreimikhailov

import com.andreimikhailov.minitheta.core._

package object brst {

  /** Fermionic variables
    */
  abstract class Thetas
  case class Th(m:Int,n:Int) extends Thetas

  case class ThOrd(m:Int, n:Int) extends Ordered[Thetas]{
    override def compare(that:Thetas):Int = (that:Thetas) match {
      case Th(p,q) => if (p==m) (q - n) else (p - m)
    }
  }

  implicit def ordTheta(x:Thetas): Ordered[Thetas] = x match {case Th(p,q) => ThOrd(p,q)}

  /** Bosonic variables
    */
  abstract class Bosons
  /** This bosonic variable is treated specially. It is reserved to be converted
    * into the unknown variable of the linear equations. The convention is:
    * (Cf->n) goes into z(n)
    * (This trick, obviously, works only for '''linear''' equations)
    */
  case class Cf() extends Bosons


  abstract class EqCoeff
  case class z(m:Int) extends EqCoeff
  def varList(start:Int,end:Int) =
    for (i <- List.range(start,end)) yield EBLS[BigInt,EqCoeff](z(i))


  case class const() extends EqCoeff {
    override def toString = "#"
  }

  /** Super-monomial
    */
  type BosonsFermions = Tuple2[Set[Thetas], Map[Bosons,Int]]
  /** Linear combination of super-monomials
    */
  type SUSY = ESCA[BigInt,Thetas,Bosons]

  implicit def ESCAtoELS(x: SUSY) :ELS[BigInt, BosonsFermions] = x match {case ESCA(y) => y}

  implicit def ERESCAtoESCA(x:ER[SUSY]):SUSY = x match { case y@ ESCA(_) => y }

  def esca(x:ELS[BigInt,BosonsFermions]) = ESCA[BigInt,Thetas,Bosons](x)
  def monoBase(t:Set[Thetas], m:Map[Bosons,Int]):EBLS[BigInt, BosonsFermions] =
    EBLS[BigInt, BosonsFermions]((t,m))
  def mono(t:Set[Thetas], m:Map[Bosons,Int]):ESCA[BigInt,Thetas,Bosons] =
    esca(monoBase(t,m))


  /** Generates unknowns z(n) in equations.
    * Checks if the Map of e contains Cf to the power n
    * If it does, erase this Cf from the map and instead inject z(n)
    * in a coefficient in front.
    */
  private def PairOfCoefAndMon (c: ER[BigInt], e: BosonsFermions)
    : Tuple2[ELS[BigInt,EqCoeff], BosonsFermions] =
      if (e._2.keys.toList.contains(Cf()))
        (c *: EBLS[BigInt,EqCoeff](z(e._2(Cf()))), ( e._1 , e._2 - Cf() ))
      else ( c *: EBLS[BigInt,EqCoeff](const()), e )

  private def collectPairs(x: SUSY)
  : List[Tuple2[ELS[BigInt,EqCoeff], BosonsFermions]] =
    x match {
        case ESCA(STimesV(s,EBLS(v))) => List(PairOfCoefAndMon(s,v))
        case ESCA(EBLS(v)) => List(PairOfCoefAndMon(UnitOfRing[BigInt](1),v))
        case ESCA(VSum(vs@ _*))       => vs.toList map {
          case STimesV(s,EBLS(v)) => PairOfCoefAndMon(s,v)
          case EBLS(v) => PairOfCoefAndMon(UnitOfRing[BigInt](1), v)
        }}
  /** Sums the coefficients of every monomial (element of BosonsFermions).
    * Returns the list of the sums of coefficients for every monomials.
    * Note that every sum has terms containing the unknowns z(n) and also terms
    * not containing z(n).
    */
  private def makeEqsFromPairs(
    xs: List[Tuple2[ELS[BigInt,EqCoeff], BosonsFermions]]
  ) : List[ELS[BigInt,EqCoeff]] = {
    val monomials = Set( (xs map (_._2)): _*).toList
    monomials map ( mnml =>
      (
        ((ZeroV[BigInt, EqCoeff](): ELS[BigInt,EqCoeff]) /:
         (xs filter (_._2 == mnml)))
        (_ + _._1)
      ))
  }
  /** Collects equations for vanishing of x,
    * i.e. vanishing of the coefficient of every monomial in x.
    * (Monomial is an element of BosonsFermions.)
    * Returns the list of equations. Every equation is a linear combination
    * of terms. There are 2 types of terms: the constants (rational numbers)
    * and rational number times z(n), where z(n) are unknowns to be solved for.
    */
  def collectEquations(x: SUSY) : List[ELS[BigInt,EqCoeff]] =
    makeEqsFromPairs(collectPairs(x))

  def unknown(n: Int)    = mono(Set(), Map(Cf().asInstanceOf[Bosons] -> n))
  def scalar(n:BigInt)   = UnitOfRing[BigInt](n)
  def ratio(m:BigInt, n:BigInt) = EmbedQ[BigInt](new Rational(m,n))
  def scalarESCA(n:BigInt) = esca(scalar(n) *: monoBase(Set(),Map()))
  def ratioESCA(m:BigInt, n:BigInt) =
    esca(ratio(m,n) *: monoBase(Set(),Map()))


  def mscalar(m:BigInt, n:BigInt, size:Int) = 
    Matrix.unit(esca(ratio(m,n) *: monoBase(Set(),Map())),
                esca(ZeroV()), 
                size) 
  val dim = 4
  val mTh = Mat[ESCA[BigInt, Thetas, Bosons]](
      for (i <- List.range(0,dim)) 
        yield for (j <- List.range(0,dim)) 
          yield
          ESCA[BigInt, Thetas, Bosons](EBLS((Set(Th(i,j).asInstanceOf[Thetas]), Map())))
      )
  val cap = Matrix.symplectic[ESCA[BigInt, Thetas, Bosons]](
      ESCA[BigInt,Thetas,Bosons](EBLS((Set(), Map()))),
      ESCA[BigInt,Thetas,Bosons](ZeroV()), dim/2)
  val cup = mscalar(-1,1,dim)*cap
  def traceless(x: Mat[ESCA[BigInt, Thetas, Bosons]])
  : Mat[ESCA[BigInt, Thetas, Bosons]] = x + (
      Matrix.unit(esca(ratio(-1, dim) *: monoBase(Set(), Map())),
                  esca(ZeroV()), 
                  dim) * 
      Matrix.unit(x.trace, ESCA[BigInt,Thetas,Bosons](ZeroV()), dim)
      )
  def cupless(x: Mat[ESCA[BigInt, Thetas, Bosons]]) 
  : Mat[ESCA[BigInt, Thetas, Bosons]] = 
    x + (mscalar(-1,dim,dim) * 
         Matrix.unit((x*cap).trace, esca(ZeroV()), dim) * 
         cup)
  def capless(x: Mat[ESCA[BigInt, Thetas, Bosons]]) 
  : Mat[ESCA[BigInt, Thetas, Bosons]] = 
    x + (mscalar(-1,dim,dim) * 
         Matrix.unit((x*cup).trace, esca(ZeroV()), dim) * 
         cap)  

  def QL(x:ESCA[BigInt,Thetas,Bosons]) :ESCA[BigInt,Thetas,Bosons] =  (
      ESCA[BigInt,Thetas,Bosons](ZeroV():ELS[BigInt,BosonsFermions]) /: 
      (for (i <- List.range(0,dim); j <- List.range(0,dim)) 
        yield (Lambda()(i)(j) * (Th(i,j) D_: x)) ) ) (_ + _)
        
  def QR(x:ESCA[BigInt,Thetas,Bosons]) :ESCA[BigInt,Thetas,Bosons] =  (
      ESCA[BigInt,Thetas,Bosons](ZeroV():ELS[BigInt,BosonsFermions]) /: 
      (for (i <- List.range(0,dim); j <- List.range(0,dim)) 
        yield (Mu()(i)(j) * (Th(i,j) D_: x)) ) ) (_ + _)
        
  def QLmat(x: Mat[ESCA[BigInt, Thetas, Bosons]]) 
  : Mat[ESCA[BigInt, Thetas, Bosons]] = Mat[ESCA[BigInt, Thetas, Bosons]](
      for (i <- List.range(0,dim)) 
        yield for (j <- List.range(0,dim)) 
          yield
          QL(x.elemss(i)(j))
      ) 
      
  def QRmat(x: Mat[ESCA[BigInt, Thetas, Bosons]]) 
  : Mat[ESCA[BigInt, Thetas, Bosons]] = Mat[ESCA[BigInt, Thetas, Bosons]](
      for (i <- List.range(0,dim)) 
        yield for (j <- List.range(0,dim)) 
          yield
          QR(x.elemss(i)(j))
      )
      
  val mLa = Mat[ESCA[BigInt,Thetas,Bosons]](Lambda())
  val mMu = Mat[ESCA[BigInt,Thetas,Bosons]](Mu())
      
  def word(wrd: String) : Mat[ESCA[BigInt, Thetas, Bosons]] = {
      def flit(c:Char) = c match {
        case '_' => cup
        case '^' => cap
        case 'l' => mLa
        case 'm' => mMu
        case 't' => mTh
      }
      def fpcn(previous: Option[Char], c:Char, next: Option[Char]) = 
        (previous,c,next) match {
          case (Some('_'), c1, _ )    => flit(c1)
          case (Some('^'), c1, _ )    => flit(c1).transpose
          case (None, c1, Some('^'))  => flit(c1)
          case (None, c1, Some('_'))  => flit(c1).transpose
          case (_, c1, _)             => flit(c1)
      }
      def main(cs: List[Char], acc:List[Char])
      : Mat[ESCA[BigInt, Thetas, Bosons]] = cs match {
        case List()  => mscalar(1,1,dim)
        case List(c) => 
          if (acc == List()) flit(c) else fpcn(Some(acc.head), c, None)
        case c::rest => 
          if (acc == List()) 
            (fpcn(None, c, Some(rest.head)) * main(rest, c::acc)) 
          else     
            (fpcn(Some(acc.head), c, Some(rest.head)) * main(rest, c::acc))
      }
      main(wrd.toList, List())
  }

  //let us verify that the Pure spinor constraints hold:
  assert(
    cupless(word("l^l")) == 
    Mat(
      for (i<-List.range(0,dim)) yield
        for (j<-List.range(0,dim)) yield
          esca(ZeroV())))
  assert(
    cupless(word("m^m")) == 
    Mat(
      for (i<-List.range(0,dim)) yield
        for (j<-List.range(0,dim)) yield
          esca(ZeroV())))
  
  def matScalarUnknown(n: Int) = Mat(for (i<-List.range(0,dim))
    yield for (j<-List.range(0,dim))
      yield if (i==j) unknown(n) else ESCA[BigInt,Thetas,Bosons](ZeroV()))
  def matScalarESCA(n: Int) = Mat(for (i<-List.range(0,dim))
    yield for (j<-List.range(0,dim))
      yield if (i==j) scalarESCA(n) else ESCA[BigInt,Thetas,Bosons](ZeroV()))
  def matScalar(x:ESCA[BigInt,Thetas,Bosons]) = Mat(for (i<-List.range(0,dim))
    yield for (j<-List.range(0,dim))
      yield if (i==j) x else ESCA[BigInt,Thetas,Bosons](ZeroV()))
  def collectEquationsFromMat(x: Mat[ESCA[BigInt, Thetas, Bosons]]) 
  : List[ELS[BigInt,EqCoeff]] =
    x.elemss.flatten flatMap (u => collectEquations(u))

  /** This is $(\theta\Gamma_m\lambda)\Gamma_m\theta$
    * and $(\theta\Gamma_m\mu)\Gamma_m\theta$
    */
  val TTL = (scalarESCA(2)*:word("t^l_t")) + word("t^t_l") + word("l^t_t")
  val TTM = (scalarESCA(2)*:word("t^m_t")) + word("t^t_m") + word("m^t_t")
  assert(QLmat(TTL) == mscalar(0,1,4))
  /** This is the AdS part of the vector class
    */
  val TLup = traceless(word("t^l_") + word("l^t_"))
  val TMup = traceless(word("t^m_") + word("m^t_"))
  /** This is the sphere part of the vector class
    */
  val TLdn = traceless(word("^t_l") + word("^l_t"))
  val TMdn = traceless(word("^t_m") + word("^m_t"))

  def insUnknownsMat(n: Int, xs: List[Mat[ESCA[BigInt, Thetas, Bosons]]]) = {
    val ys = for (i <- List.range(0, xs.length)) yield unknown(n + i) *: xs(i)
    (mscalar(0, 1, dim) /: ys)(_ + _)
  }

 }