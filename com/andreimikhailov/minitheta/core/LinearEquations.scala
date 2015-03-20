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

import annotation.tailrec

object EquationHelpers {
  def substitute[S,T] (m: Map[EBLS[S,T], ELS[S,T]], 
                       expression: ELS[S,T], 
                       whichOnes: List[EBLS[S,T]] ) : ELS[S,T] =
      ((ZeroV[S,T]():ELS[S,T]) /: (expression.monomials map { 
        // case ZeroV()      => ZeroV[S,T]()
        case x@ EBLS(_)   => 
          if ((m.keys.toList.contains(x)) && (whichOnes.contains(x))) m(x) 
          else x
        case STimesV(s,x) => 
          if ((m.keys.toList.contains(x)) && (whichOnes.contains(x))) s *: m(x) 
          else STimesV[S,T](s,x)
                      }))    (_ + _) 
}
/** This solves a single linear equation
  * A linear equation is a sum of monomials, some of which do not depend on x (constant terms)
  * and others are proportional to x. Let us look at the first monomial. If it is proportional
  * to x, then we resolve for that x. If it does not depend on x, move it to the accumulator
  * and look at the next term.
  */
case class Equation[S,T](shouldBeZero: ELS[S,T], variables: List[EBLS[S,T]]){ 
  @tailrec private def solve_(eqn: Equation[S,T], 
                              acc:List[ELS[S,T] with IsMonomial]) 
  : Option[Map[EBLS[S,T], ELS[S,T]]] = 
    {    
    eqn.shouldBeZero match { 
      case ZeroV()        => Some(Map())
      case _              => eqn.shouldBeZero.monomials match {
        case List()               => 
          if (acc == List()) Some(Map()) 
          else None
        case List(el@ EBLS(_))    => if (variables.contains(el)) 
          Some(Map(el -> UnitOfRing[S](-1) *: (
                  ((ZeroV():ELS[S,T]) /: acc) (_ + _)
                ))) 
          else None
        case List(STimesV(x,v))   => 
          if (variables.contains(v))  
            Some(Map(v  -> 
                     (UnitOfRing[S](-1) * x.inv) *: (
                      ((ZeroV():ELS[S,T]) /: acc) (_ + _)
                    ))) 
          else None
        case (el@ EBLS(_))::rest  => 
          if (variables.contains(el)) 
            Some(Map(el -> (((ZeroV():ELS[S,T]) /: (rest:::acc)) 
                            ((u,v) => u +  (UnitOfRing[S](-1) *: v))
                )))
          else solve_( 
            Equation[S,T]((((ZeroV[S,T]():ELS[S,T]) /: rest) (_ + _)), 
                          variables), 
            el::acc
          )
        case (STimesV(x,v))::rest => 
          if (variables.contains(v))  
            Some(Map(v  -> (((ZeroV():ELS[S,T]) /: (rest:::acc)) 
                            ((u,v) => u + ((UnitOfRing[S](-1) * x.inv) *: v))
                )))
            else solve_( 
              Equation[S,T]((((ZeroV[S,T]():ELS[S,T]) /: rest) (_ + _) ), 
                            variables), 
              STimesV[S,T](x,v)::acc
            )
      }
    }
  }
  /** This method returns either None or Some map of either 0 or 1 keys
    * None is returned when there is no solution (incompatible system)
    * Empty Map (''i.e.'' Map()) means the opposite: no constraints
    */
  def solve = {
    val soln = solve_(this, List()) 
    require (soln match {
        case None => true case Some(mp) => (mp.keys.toList.length < 2)
      } )
    soln
  }
  def substitute ( m: Map[EBLS[S,T], ELS[S,T]] ) : ELS[S,T] = 
    EquationHelpers.substitute[S,T](m, shouldBeZero, variables) 
}
case class Equations[S,T](eqs: List[ELS[S,T]], variables: List[EBLS[S,T]]) { 
  type RULES = Map[EBLS[S,T],ELS[S,T]]
  /** This method solves a system of equations keeping track of the order
    * in which the unknowns were solved. 
    * (The second element List in Tuple2 is to keep track of the order
    * of resolved variables.)
    */
  private def solve1   //TODO: make it tailrec!
  : (List[ELS[S,T]] => Tuple2[ Option[RULES], List[EBLS[S,T]] ]) = {
    case List()        => (Some(Map()), List())
    case List(ZeroV()) => (Some(Map()), List())
    case List(eq)      => {
      val res = Equation[S,T](eq, variables).solve
      res match {
        case None => (None, List()) 
        case Some(r) => {require(r.keys.toList.length < 2) ;
                         (Some(r),r.keys.toList)}
        } 
    }
    case eq::rest  => {
      val m1: Option[RULES] = Equation[S,T](eq, variables).solve
      m1 match {
        case None => (None,List())
        case Some(rl) => {
          // Note that rl is either Map() or Map(a->b), i.e. rl.keys.toList.length < 2
          assert(rl.keys.toList.length < 2)
          val restSub = 
            rest map (a => Equation[S,T](a, variables).substitute(rl))
          val solnRest:Tuple2[Option[RULES], List[EBLS[S,T]]] = solve1(restSub)
          solnRest match {
            case (None,_)             => (None, List())
            case (Some(mrest), lrest) => {
              require(rl.keys.toList match {
                  case List()   => true 
                  case List(r1) => !(mrest.keys.toList.contains(r1))})
              (Some(rl ++ mrest), lrest:::(rl.keys.toList)) 
            }}}}}} 
  @tailrec private def reduceRules(
    ruleMap: Map[EBLS[S,T], ELS[S,T]], 
    basisEls: List[EBLS[S,T]],
    acc: Map[EBLS[S,T], ELS[S,T]]
  ): Map[EBLS[S,T], ELS[S,T]] = { 
    require (ruleMap.keySet == Set(basisEls: _*)) // basisEls is just to keep track of the ordering of the rules
    require (Set(basisEls: _*).count(_ => true) == basisEls.length) // basisEls should be all different
    basisEls match {
      case List()    => ruleMap ++ acc
      case List(b)   => ruleMap ++ acc
      case b::bs     => reduceRules(
          Map( 
            (for (x <- bs) 
              yield (
                  x -> EquationHelpers.substitute(
                    Map(b -> ruleMap(b)), 
                    ruleMap(x), 
                    List(b)
                  ) )
            ) : _*  
          ) ,   
          bs,
          Map(b-> ruleMap(b)) ++ acc
        )
    }
  }  
  /** Solves the system of equations:
    */
  def solve: Option[RULES] = solve1(this.eqs) match {
    case (None,_)                  => None
    case (Some(ruleMap), basisEls) => {
      Some(reduceRules(ruleMap, basisEls, Map()))
    }}} 
