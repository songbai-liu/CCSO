//  MaF13.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.problems.MaF;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

import java.util.Vector;

/**
 * Class representing problem MaF13
 */
public class MaF13 extends Problem {

	/**
	 * Creates a default MaF13 problem (5 variables and M objectives)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF13(String solutionType) throws ClassNotFoundException {
		this(solutionType, 5, 4);
	} // LZ09_F6

	/**
	 * Creates a LZ09_F6 problem instance
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF13(String solutionType, Integer numberOfVariables,
			Integer numberOfObjectives) {
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = numberOfObjectives;
		numberOfConstraints_ = 0;
		problemName_ = "MaF13";

		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		for (int var = 0; var < 2; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for
		
		for (int var = 2; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -2.0;
			upperLimit_[var] = 2.0;
		} // for

		if (solutionType.compareTo("BinaryReal") == 0)
			solutionType_ = new BinaryRealSolutionType(this);
		else if (solutionType.compareTo("Real") == 0)
			solutionType_ = new RealSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType
					+ " invalid");
			System.exit(-1);
		}
	} // MaF13

	/**
	 * Evaluates a solution
	 * 
	 * @param solution
	 *            The solution to evaluate
	 * @throws JMException
	 */
	public void evaluate(Solution solution) throws JMException {
		Variable[] gen = solution.getDecisionVariables();

		Vector<Double> x = new Vector<Double>(numberOfVariables_);
		Vector<Double> y = new Vector<Double>(numberOfObjectives_);

		for (int i = 0; i < numberOfVariables_; i++) {
			x.addElement(gen[i].getValue());
			y.addElement(0.0);
		} // for
		
		for (int i = 0; i < numberOfObjectives_; i++) {
			y.addElement(0.0);
		} // for


		Vector<Double> aa = new Vector();
		Vector<Double> bb = new Vector();
		Vector<Double> cc = new Vector();
		Vector<Double> dd = new Vector();
		
		double y_value;
		
		for (int n = 2; n < numberOfVariables_; n++) {
			y_value = x.elementAt(n) - 2*x.elementAt(1)*Math.sin(2*x.elementAt(0)*Math.PI+n*Math.PI/numberOfVariables_);
			if (n % 3 == 1)
				aa.addElement(y_value);
			else if (n % 3 == 2)
				bb.addElement(y_value);
			else if(n % 3 == 0){
				cc.addElement(y_value);
			}
			if(n>=3){
				dd.addElement(y_value);
			}
		}
		double a=0.0;
		double b=0.0;
		double c=0.0;
		double d=0.0;
		for(int i=0;i<aa.size();i++){
			a += Math.pow(aa.get(i), 2);
		}
		for(int i=0;i<bb.size();i++){
			b += Math.pow(bb.get(i), 2);
		}
		for(int i=0;i<cc.size();i++){
			c += Math.pow(cc.get(i), 2);
		}
		for(int i=0;i<dd.size();i++){
			d += Math.pow(dd.get(i), 2);
		}
		a *= 2.0/aa.size();
		b *= 2.0/bb.size();
		c *= 2.0/cc.size();
		d *= 2.0/dd.size();
		
		double a_value = Math.sin(0.5*Math.PI*x.elementAt(0)) + a;
		double b_value = Math.cos(0.5*Math.PI*x.elementAt(0))*Math.sin(0.5*Math.PI*x.elementAt(1)) + b;
		double c_value = Math.cos(0.5*Math.PI*x.elementAt(0))*Math.cos(0.5*Math.PI*x.elementAt(1)) + c;
		double d_value = Math.pow(a_value, 2)+Math.pow(b_value, 10)+Math.pow(c_value, 10)+d;
		
		y.set(0, a_value);
		y.set(1, b_value);
		y.set(2, c_value);
		
		for (int i = 3; i < numberOfObjectives_; i++)
			y.set(i, d_value);
		
		for (int i = 0; i < numberOfObjectives_; i++)
			solution.setObjective(i, y.get(i));
		
		
	} // evaluate
} //MaF13

