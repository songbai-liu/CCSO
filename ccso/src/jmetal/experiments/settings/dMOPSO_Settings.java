//  dMOEAD_Settings.java
//
//  Authors:
//       Jorge Rodriguez
//       Antonio J. Nebro <antonioo@lcc.uma.es>
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

package jmetal.experiments.settings;

import jmetal.core.Algorithm;
import jmetal.experiments.Settings;
import jmetal.metaheuristics.dmopso.dMOPSO;
import jmetal.problems.ProblemFactory;
import jmetal.util.JMException;

/**
 * Settings class of algorithm dMOPSO
 */
public class dMOPSO_Settings extends Settings {
	private String dataDirectory_;
	private int swarmSize_;
	private int maxIterations_;
	private int maxAge_;
	private int funcIntType_;
	private String functionType_;

	/**
	 * Constructor
	 */
	public dMOPSO_Settings(String problem) {
		super(problem);

		Object[] problemParams = { "Real" };
		try {
			problem_ = (new ProblemFactory()).getProblem(problemName_,
					problemParams);
		} catch (JMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		// Default experiments.settings
		swarmSize_ = 100;
		maxIterations_ = 250;
		maxAge_ = 2;
		functionType_ = "";
		funcIntType_ = 0;

		// Directory with the files containing the weight vectors used in
		// Q. Zhang, W. Liu, and H Li, The Performance of a New Version of
		// MOEA/D
		// on CEC09 Unconstrained MOP Test Instances Working Report CES-491,
		// School
		// of CS & EE, University of Essex, 02/2009.
		// http://dces.essex.ac.uk/staff/qzhang/MOEAcompetition/CEC09final/code/ZhangMOEADcode/moead0305.rar

		dataDirectory_ = "/Users/antelverde/Softw/pruebas/data/MOEAD_parameters/Weight";
		// dataDirectory_ = "/home/jorgero/moead/weight";

	} // MOEAD_Settings

	/**
	 * Configure the algorithm with the specified parameter experiments.settings
	 * 
	 * @return an algorithm object
	 * @throws jmetal.util.JMException
	 */
	public Algorithm configure() throws JMException {
		Algorithm algorithm;

		// Creating the problem
		algorithm = new dMOPSO(problem_);

		switch (funcIntType_) {
		case 0:
			functionType_ = "_PBI";
			break;
		case 1:
			functionType_ = "_TCHE";
			break;
		case 2:
			functionType_ = "_AGG";
			break;
		}

		// Algorithm parameters
		algorithm.setInputParameter("swarmSize", swarmSize_);
		algorithm.setInputParameter("maxIterations", maxIterations_);
		algorithm.setInputParameter("maxAge", maxAge_);
		algorithm.setInputParameter("functionType", functionType_);
		algorithm.setInputParameter("dataDirectory", dataDirectory_);

		return algorithm;
	} // configure
} // MOEAD_Settings
