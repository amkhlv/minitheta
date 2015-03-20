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

package com.andreimikhailov.brst

import com.andreimikhailov.minitheta.core._
import com.andreimikhailov.minitheta.examples._

case class Alpha() extends Bosons
case class Beta() extends Bosons
case class Gamma() extends Bosons

/** Lambda and Mu are the canonical form of a pair of pure spinors
  * in generic position
  */
object Lambda {
  def apply() : List[List[ESCA[BigInt,Thetas,Bosons]]] = List(
      List('a', 0 , 0 , 0 ), 
      List( 0 ,'b', 0 , 0 ),
      List( 0 , 0 ,'c', 0 ),
      List( 0 , 0 , 0 ,'d')
      ) map {x => x map {
          case 0   => esca(ZeroV())
          case 1   => mono(Set(),Map())
          case 'a' => mono(Set(), Map(Alpha() -> 1))
          case 'b' => mono(Set(), Map(Beta()  -> 1, Gamma() -> 1))
          case 'c' => mono(Set(), Map(Gamma() -> 1))
          case 'd' => mono(Set(), Map(Alpha() -> 1, Beta()  -> 1))
          }}}

/** Lambda and Mu are the canonical form of a pair of pure spinors
  * in generic position
  */
object Mu {
  def apply() : List[List[ESCA[BigInt,Thetas,Bosons]]] = List(
      List(1,  0,  0,  0), 
      List(0,  1,  0,  0),
      List(0,  0,  1,  0),
      List(0,  0,  0,  1)
      ) map {x => x map {
          case 0 => esca(ZeroV())
          case z => esca(scalar(z) *: monoBase(Set(),Map()))
          }}}
