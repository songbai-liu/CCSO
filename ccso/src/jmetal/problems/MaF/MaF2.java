//  MaF2.java
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
 * Class representing problem MaF2
 */
public class MaF2 extends Problem {

	/**
	 * Creates a default MaF2 problem (12 variables and 3 objectives)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF2(String solutionType) throws ClassNotFoundException {
		this(solutionType, 12, 3);
	} // DTLZ2

	/**
	 * Creates a new instance of MaF2
	 * 
	 * @param numberOfVariables
	 *            Number of variables
	 * @param numberOfObjectives
	 *            Number of objective functions
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF2(String solutionType, Integer numberOfVariables,
			Integer numberOfObjectives) {
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = numberOfObjectives;
		numberOfConstraints_ = 0;
		problemName_ = "MaF2";

		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];
		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
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
	} // MaF2

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
		double[] g = new double[numberOfObjectives_];
		double[] theta = new double[numberOfObjectives_-1];
		int k = numberOfVariables_ - numberOfObjectives_ + 1;

		for (int i = 0; i < numberOfVariables_; i++)
			x[i] = gen[i].getValue();
		
		for(int i = 1; i <= numberOfObjectives_-1; i++){
			int down = numberOfObjectives_+(i-1)*(int)Math.floor((numberOfVariables_-numberOfObjectives_+1)/numberOfObjectives_)-1;
			int up = (numberOfObjectives_)+(i)*(int)Math.floor((numberOfVariables_-numberOfObjectives_+1)/numberOfObjectives_)-2;
			//System.out.println("down = "+ down + ", up = "+up);
			double value = 0.0;
			for(int j=down;j<=up;j++){
				value += Math.pow(((x[j]*0.5 + 0.25)-0.5), 2.0);
			}
			g[i-1] = value;
		}
		int d = numberOfObjectives_ + (numberOfObjectives_-1)*(int)Math.floor((numberOfVariables_-numberOfObjectives_+1)/numberOfObjectives_) - 1;
		//System.out.println("d = "+d);
		double v = 0.0;
		for(int i=d;i<numberOfVariables_;i++){
			v += Math.pow(((x[i]*0.5 + 0.25)-0.5), 2.0);
		}
		g[numberOfObjectives_-1] = v;

		for (int i = 0; i < numberOfObjectives_; i++)
			f[i] = 1.0 + g[i];
		
		for (int i = 0; i < numberOfObjectives_-1; i++) {
			theta[i] = (0.5*(Math.PI))*(x[i]*0.5 + 0.25);
		}

		for (int i = 0; i < numberOfObjectives_; i++) {
			for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
				f[i] *= Math.cos(theta[j]);
			if (i != 0) {
				int aux = numberOfObjectives_ - (i + 1);
				f[i] *= Math.sin(theta[aux]);
			} // if
		} // for

		for (int i = 0; i < numberOfObjectives_; i++)
			solution.setObjective(i, f[i]);
	}
} // evaluate
