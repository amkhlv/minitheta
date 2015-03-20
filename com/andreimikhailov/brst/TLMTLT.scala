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
import  com.andreimikhailov.minitheta.core._
import  com.andreimikhailov.minitheta.examples._

/** This verifies the following equation from [[http://arxiv.org/abs/1105.2231 arXiv:1105.2231]] :
  * <br>
  * <img src="http://andreimikhailov.com/media/minitheta-scaladoc/Eq-181.png" alt="Eq. (181)">
  * <br>
  * This equation holds in cohomology of QL, ''i.e.'' up to QL-exact terms
  */

object TLMTLT extends App {

  /** This is
    *   <img src="http://andreimikhailov.com/media/minitheta-scaladoc/TLTTLT.png" alt="Lambda2Theta4">
    */
  val TLTTLT  = cupless(word("t^l_t^t_l^t"))

  /** This is
    *  <img src="http://andreimikhailov.com/media/minitheta-scaladoc/QrTLTTLT.png" alt="QrLambda2Theta4">
    */
  val QrTLTTLT =
    QRmat(cupless(word("t^l_t^t_l^t")))

  /** This is
    * <img src="http://andreimikhailov.com/media/minitheta-scaladoc/MTLT_tr-times-LT.png" alt="MTLT_tr-times-LT">
    */

  val MTLTnrmLT =
    (word("m^t_l^t_").trace) *: cupless(word("t^l") + word("l^t"))

  /** This is
    * <img src="http://andreimikhailov.com/media/minitheta-scaladoc/MLTTLT.png" alt="MLTTLT">
    */
  val RHS =
    cupless(word("m^l_t^t_l^t") + (scalarESCA(-1) *: word("t^l_t^t_l^m"))) +
  ((ratioESCA(-1,4) * word("l^t_").trace) *: cupless(word("m^t_l^t"))) +
  ((ratioESCA(1,4)  * word("t^l_").trace) *: cupless(word("t^l_t^m")))
      
  /** This is the ansatz for the gauge transformation, ''i.e.'', the equation which we are verifying holds
    * up to QL(RHS)
    */
  val MLT4   = 
    (unknown(100) *: word("m^l_t^t_t^t")) +
    (unknown(101) *: word("m^t_l^t_t^t")) +
    (unknown(102) *: word("m^t_t^l_t^t")) +
    (unknown(103) *: word("m^t_t^t_l^t")) +
    (unknown(104) *: word("m^t_t^t_t^l")) +
    (unknown(105) *: word("l^m_t^t_t^t")) +
    (unknown(106) *: word("t^m_l^t_t^t")) +
    (unknown(107) *: word("t^m_t^l_t^t")) +
    (unknown(108) *: word("t^m_t^t_l^t")) +
    (unknown(109) *: word("t^m_t^t_t^l")) +
    (unknown(110) *: word("l^t_m^t_t^t")) +
    (unknown(111) *: word("t^l_m^t_t^t")) +
    (unknown(112) *: word("t^t_m^l_t^t")) +
    (unknown(113) *: word("t^t_m^t_l^t")) +
    (unknown(114) *: word("t^t_m^t_t^l")) +
    (unknown(115) *: word("l^t_t^m_t^t")) +
    (unknown(116) *: word("t^l_t^m_t^t")) +
    (unknown(117) *: word("t^t_l^m_t^t")) +
    (unknown(118) *: word("t^t_t^m_l^t")) +
    (unknown(119) *: word("t^t_t^m_t^l")) +
    (unknown(120) *: word("l^t_t^t_m^t")) +
    (unknown(121) *: word("t^l_t^t_m^t")) +
    (unknown(122) *: word("t^t_l^t_m^t")) +
    (unknown(123) *: word("t^t_t^l_m^t")) +
    (unknown(124) *: word("t^t_t^t_m^l")) +
    (unknown(125) *: word("l^t_t^t_t^m")) +
    (unknown(126) *: word("t^l_t^t_t^m")) +
    (unknown(127) *: word("t^t_l^t_t^m")) +
    (unknown(128) *: word("t^t_t^l_t^m")) +
    (unknown(129) *: word("t^t_t^t_l^m")) +
    ((unknown(130) * word("t^t_t^t_").trace) *: word("l^m")) +
    ((unknown(131) * word("t^t_t^t_").trace) *: word("m^l")) +
    ((unknown(132) * word("l^t_t^t_").trace) *: word("m^t")) +
    ((unknown(133) * word("l^t_t^t_").trace) *: word("t^m")) +
    ((unknown(134) * word("m^t_t^t_").trace) *: word("l^t")) +
    ((unknown(135) * word("m^t_t^t_").trace) *: word("t^l")) +
    ((unknown(136) * word("l^m_t^t_").trace) *: word("t^t")) +
    ((unknown(137) * word("l^t_m^t_").trace) *: word("t^t")) +
    ((unknown(138) * word("l^t_t^m_").trace) *: word("t^t")) +
    ((unknown(139) * (word("l^t_").trace)*(word("m^t_").trace)) *: 
     word("t^t")) +
    ((unknown(140) * (word("l^m_").trace)) *: word("t^t_t^t")) +
    ((unknown(141) * (word("l^t_").trace)) *: word("m^t_t^t")) +
    ((unknown(142) * (word("l^t_").trace)) *: word("t^m_t^t")) +
    ((unknown(143) * (word("l^t_").trace)) *: word("t^t_m^t")) +
    ((unknown(144) * (word("l^t_").trace)) *: word("t^t_t^m")) +
    ((unknown(145) * (word("m^t_").trace)) *: word("l^t_t^t")) +
    ((unknown(146) * (word("m^t_").trace)) *: word("t^l_t^t")) +
    ((unknown(147) * (word("m^t_").trace)) *: word("t^t_l^t")) +
    ((unknown(148) * (word("m^t_").trace)) *: word("t^t_t^l")) +
    ((unknown(150) * word("l^m_t^t_t^t_").trace) *: cup) +
    ((unknown(151) * word("l^t_m^t_t^t_").trace) *: cup) +
    ((unknown(152) * word("l^t_t^m_t^t_").trace) *: cup) +
    ((unknown(153) * word("l^t_t^t_m^t_").trace) *: cup) +
    ((unknown(154) * word("l^t_t^t_t^m_").trace) *: cup) +
    ((unknown(155) * word("l^m_").trace * word("t^t_t^t_").trace) *: cup) +
    ((unknown(156) * word("l^t_").trace * word("m^t_t^t_").trace) *: cup) +
    ((unknown(157) * word("m^t_").trace * word("l^t_t^t_").trace) *: cup)

  /** Here z(1) is defined as follows:
    * <br>
    *   <img src="http://andreimikhailov.com/media/minitheta-scaladoc/QrTLTTLT.png" alt="QrLambda2Theta4">
    * <br>
    * plus z(1) times
    * <br>
    *   <img src="http://andreimikhailov.com/media/minitheta-scaladoc/MLTTLT.png" alt="MLTTLT">
    * <br>
    * equals QL of something
    */
  val zmap1 =  
    Equations(collectEquationsFromMat(
      QRmat(TLTTLT) + (unknown(1)*:RHS) + QLmat(MLT4)
      ), varList(0,300)
    ).solve
  

  /** Here z(2) is defined as follows:
    * <br>
    *   <img src="http://andreimikhailov.com/media/minitheta-scaladoc/MTLT_tr-times-LT.png" alt="MTLT_tr-times-LT">
    * <br>
    *   plus z(2) times
    * <br>
    *   <img src="http://andreimikhailov.com/media/minitheta-scaladoc/MLTTLT.png" alt="MLTTLT">
    * <br>
    * equals QL of something. From [[http://arxiv.org/abs/1105.2231 arXiv:1105.2231]] :
    * <br>
    *   <img src="http://andreimikhailov.com/media/minitheta-scaladoc/ActingOnTheDoubleTrace.png" alt="ActingOnTheDoubleTrace">
    */
  val zmap2 =
    Equations(collectEquationsFromMat(
      MTLTnrmLT + (unknown(2)*:RHS) + QLmat(MLT4)
      ), varList(0,300)
    ).solve    


  val zm1 = zmap1 match {
    case Some(x) => x
    case None => throw  new RuntimeException("NOTHING FOUND")}
  println("zmap1 is:")
  for (k <- zm1.keys) println(k + "--->  " + zm1(k))
  // unknown(1)--->  -8/5
  println("zmap2 is:")
  val zm2 = zmap2 match {
    case Some(x) => x 
    case None => throw  new RuntimeException("NOTHING FOUND")}
  for (k <- zm2.keys) println(k + "--->  " + zm2(k))
  // unknown(2)---> -4/5
}