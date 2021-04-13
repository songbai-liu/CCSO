//  Utils.java
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

package jmetal.metaheuristics.ccso;

import jmetal.core.Solution;

/**
 * Utilities methods to used by MCCSO
 */
public class Utils {

	public static double distVector(Solution vector1, Solution vector2) {
		int dim = vector1.getNumberOfObjectives();
		double sum = 0;
		for (int n = 0; n < dim; n++) {
			sum += (vector1.getNormalizedObjective(n) - vector2.getNormalizedObjective(n)) *
					(vector1.getNormalizedObjective(n) - vector2.getNormalizedObjective(n));
		}
		return Math.sqrt(sum);
	} // distVector
	
	public static double distVector1(Solution vector1, Solution vector2) {
		int dim = vector1.getNumberOfObjectives();
		double sum = 0;
		for (int n = 0; n < dim; n++) {
			sum += (vector1.getIthTranslatedObjective(n) - vector2.getIthTranslatedObjective(n)) *
					(vector1.getIthTranslatedObjective(n) - vector2.getIthTranslatedObjective(n));
		}
		return Math.sqrt(sum);
	} // distVector
	
	public static double distVector2(double[] vector1, Solution vector2) {
		int dim = vector2.getNumberOfObjectives();
		double sum = 0;
		for (int n = 0; n < dim; n++) {
			sum += (vector1[n] - vector2.getIthTranslatedObjective(n)) *(vector1[n] - vector2.getIthTranslatedObjective(n));
		}
		return Math.sqrt(sum);
	} // distVector
	
	public static double distVector3(Solution vector1, Solution vector2) {
		int dim = vector1.getNumberOfObjectives();
		double norm = 1.0;
		double dist1 = 0.0;
		double dist2 = 0.0;
		for (int n = 0; n < dim; n++) {
			norm += (vector1.getNormalizedObjective(n)-1.0) * (vector2.getNormalizedObjective(n)-1.0);
			dist1 += (vector1.getNormalizedObjective(n)-1.0) * (vector1.getNormalizedObjective(n)-1.0);
			dist2 += (vector2.getNormalizedObjective(n)-1.0) * (vector2.getNormalizedObjective(n)-1.0);		
		}
		dist1 = Math.sqrt(dist1);
		dist2 = Math.sqrt(dist2);
		double value = Math.abs(norm/(dist1*dist2));
		if(value > 1) {
			value = 1;
		}
		if(value < 0) {
			value = 0;
		}
		value = Math.acos(value);
		return value;
	} // distVector

	public static void minFastSort(double x[], int idx[], int n, int m) {
		for (int i = 0; i < m; i++) {
			for (int j = i + 1; j < n; j++) {
				if (x[i] > x[j]) {
					double temp = x[i];
					x[i] = x[j];
					x[j] = temp;
					int id = idx[i];
					idx[i] = idx[j];
					idx[j] = id;
				} // if
			}
		} // for

	} // minFastSort

	public static void randomPermutation(int[] perm, int size) {
		int[] index = new int[size];
		boolean[] flag = new boolean[size];

		for (int n = 0; n < size; n++) {
			index[n] = n;
			flag[n] = true;
		}

		int num = 0;
		while (num < size) {
			int start = jmetal.util.PseudoRandom.randInt(0, size - 1);
			// int start = int(size*nd_uni(&rnd_uni_init));
			while (true) {
				if (flag[start]) {
					perm[num] = index[start];
					flag[start] = false;
					num++;
					break;
				}
				if (start == (size - 1)) {
					start = 0;
				} else {
					start++;
				}
			}
		} // while
	} // randomPermutation
	static void QuickSort(double[] array, int[] idx, int from, int to) {
		if (from < to) {
			double temp = array[to];
			int tempIdx = idx[to];
			int i = from - 1;
			for (int j = from; j < to; j++) {
				if (array[j] <= temp) {
					i++;
					double tempValue = array[j];
					array[j] = array[i];
					array[i] = tempValue;
					int tempIndex = idx[j];
					idx[j] = idx[i];
					idx[i] = tempIndex;
				}
			}
			array[to] = array[i + 1];
			array[i + 1] = temp;
			idx[to] = idx[i + 1];
			idx[i + 1] = tempIdx;
			QuickSort(array, idx, from, i);
			QuickSort(array, idx, i + 1, to);
		}
	}
}
