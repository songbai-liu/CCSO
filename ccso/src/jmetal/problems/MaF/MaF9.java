//  MaF9.java
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.problems.MaF;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem MaF8
 */
public class MaF9 extends Problem {

	/**
	 * Creates a default MaF8 problem instance (2 variables and 3 objectives)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF9(String solutionType) throws ClassNotFoundException {
		this(solutionType, 2, 3);
	} // DTLZ7

	/**
	 * Creates a new MaF8 problem instance
	 * 
	 * @param numberOfVariables
	 *            Number of variables
	 * @param numberOfObjectives
	 *            Number of objective functions
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF9(String solutionType, Integer numberOfVariables,
			Integer numberOfObjectives) {
		numberOfVariables_ = numberOfVariables;
		numberOfVariables_ = 2;
		numberOfObjectives_ = numberOfObjectives;
		numberOfConstraints_ = 0;
		problemName_ = "MaF9";

		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		double[][] vertexes = new double[numberOfObjectives_][2];
		vertexes[0][0] = 0.0;
		vertexes[0][1] = 1.0;
		
		double arc = 2*Math.PI/numberOfObjectives_;

		for (int i = 1; i < numberOfObjectives_; i++){
			vertexes[i][0] = vertexes[0][0] - Math.sin(arc*i);
			vertexes[i][1] = vertexes[0][1] - 1.0 + Math.cos(arc*i);
		}
		
		for (int var = 0; var < numberOfVariables_; var++) {
			//lowerLimit_[var] = -10000.0;
			//upperLimit_[var] = 10000.0;
			
			lowerLimit_[var] = -10000.0;
			upperLimit_[var] = 10000.0;
		}

		if (solutionType.compareTo("BinaryReal") == 0)
			solutionType_ = new BinaryRealSolutionType(this);
		else if (solutionType.compareTo("Real") == 0)
			solutionType_ = new RealSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType
					+ " invalid");
			System.exit(-1);
		}
	}

	/**
	 * Evaluates a solution
	 * 
	 * @param solution
	 *            The solution to evaluate
	 * @throws JMException
	 */
	public void evaluate(Solution solution) throws JMException {
		Variable[] gen = solution.getDecisionVariables();

		double[] x = new double[numberOfVariables_];
		double[] f = new double[numberOfObjectives_];

		for (int i = 0; i < numberOfVariables_; i++)
			x[i] = gen[i].getValue();
		
		double[][] point = new double[numberOfObjectives_][numberOfVariables_];
		point[0][0] = 0.0;
		point[0][1] = 1.0;
		
		double arc = 2*Math.PI/numberOfObjectives_;

		// Calculate g
		for (int i = 1; i < numberOfObjectives_; i++){
			point[i][0] = point[0][0] - Math.sin(arc*i);
			point[i][1] = point[0][1] - 1.0 + Math.cos(arc*i);
		}
        double[] k = new double[numberOfObjectives_];
		// ->Calculate fM
		for (int i = 0; i < numberOfObjectives_-1; i++){
			k[i] = (point[i+1][1]-point[i][1])/(point[i+1][0]-point[i][0]);
			f[i] = Math.abs(x[1]-k[i]*x[0]+k[i]*point[i][0]-point[i][1])/Math.sqrt(1+Math.pow(k[i], 2));
		}
		k[numberOfObjectives_-1]=(point[numberOfObjectives_-1][1]-point[0][1])/(point[numberOfObjectives_-1][0]-point[0][0]);
		f[numberOfObjectives_-1] = Math.abs(x[1]-k[numberOfObjectives_-1]*x[0]+k[numberOfObjectives_-1]*point[0][0]-point[0][1])
				/Math.sqrt(1+Math.pow(k[numberOfObjectives_-1], 2));
		// -> Setting up the value of the objetives
		for (int i = 0; i < numberOfObjectives_; i++)
			solution.setObjective(i, f[i]);
		// <-
	} // evaluate
}