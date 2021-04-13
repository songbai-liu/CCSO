package jmetal.metaheuristics.ccso;

import java.util.ArrayList;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.FitnessComparator;

public class ClusteringBasedRanking {
	/**
	 * The <code>SolutionSet</code> to rank
	 */
	private SolutionSet solutionSet_;

	/**
	 * An array containing all the fronts found during the search
	 */
	private List<SolutionSet> ranking_;
	private int clusterSize_;
	private int objNumber_;
	/**
	 * Constructor.
	 * 
	 * @param solutionSet
	 *            The <code>SolutionSet</code> to be ranked.
	 */
	public ClusteringBasedRanking(SolutionSet solutionSet, int clusterSize) {
		solutionSet_ = solutionSet;
		clusterSize_ = clusterSize;
		objNumber_ = solutionSet.get(0).getNumberOfObjectives();
		ranking_ = new ArrayList<SolutionSet>();
		ranking();
	}
	
	public void ranking() {
		SolutionSet[] sets = new MostSimilarBasedSampling(solutionSet_,objNumber_).getIdealPointOrientedPopulation(clusterSize_);
		assignFitness(sets);
		SolutionSet[] clusters = clustering(sets);
		int maxLen = Integer.MIN_VALUE;
		for (int i = 0; i < clusterSize_; i++) {
			if (clusters[i].size() > maxLen)
				maxLen = clusters[i].size();
			clusters[i].sort(new FitnessComparator());
		}
		for (int i = 0; i < maxLen; i++) {
			SolutionSet set = new SolutionSet();
			for (int j = 0; j < clusterSize_; j++) {
				if (clusters[j].size() > i) {
					clusters[j].get(i).setRank(i);
					set.add(clusters[j].get(i));
				}
			}
			ranking_.add(set);
		}
	}
	
	public SolutionSet[] clustering(SolutionSet[] sets) {
		if(sets[0].size() != clusterSize_) {
			System.out.println("the clustering size didn't equal to reference vectors");
		}
		SolutionSet[] clusters = new SolutionSet[clusterSize_];
		for(int i=0;i<clusterSize_;i++) {
			clusters[i] = new SolutionSet();
			clusters[i].add(sets[0].get(i));
			sets[0].get(i).setClusterID(i);
		}
		for(int i=0; i<sets[1].size();i++){
			Solution s1 = sets[1].get(i);
			double minDis = computeDistance(s1,sets[0].get(0));
			int minIndex = 0;
			for(int j=1; j<sets[0].size();j++){
				Solution s2 = sets[0].get(j);
				double dist = computeDistance(s1,s2);
				if(dist<minDis){
				   minDis = dist;
				   minIndex = j;
				}
			}
			clusters[minIndex].add(s1);
			s1.setClusterID(minIndex);
		}
		return clusters;
	}
	
	public double computeDistance(Solution so1, Solution so2){
		double dis = 0.0;
		double innerProduc = 0.0;
		for(int i=0; i<objNumber_; i++){
			innerProduc += Math.pow(so1.getIthTranslatedObjective(i)-so2.getIthTranslatedObjective(i), 2);
		}
		dis = Math.sqrt(innerProduc);
		return dis;
	}
	
	public void assignFitness(SolutionSet[] sets) {
		double[] lamada = new double[objNumber_];
		for(int i=0;i<objNumber_;i++) {
			lamada[i] = 0.0;
		}
		for(int i=0;i<objNumber_;i++) {
			for(int j=0;j<sets[0].size();j++) {
				lamada[i] += sets[0].get(j).getIthTranslatedObjective(i);
			}
			lamada[i] = lamada[i]/sets[0].size();
		}
		
		for(int i=0;i<sets[0].size();i++) {
			double sum = 0.0;
			for(int j=0;j<objNumber_;j++) {
				sum += sets[0].get(i).getNormalizedObjective(j)*lamada[j];
			}
			//sets[0].get(i).setFitness(sum);
			sets[0].get(i).setFitness(sets[0].get(i).getSumValue());
		}
		
		for(int i=0;i<sets[1].size();i++) {
			double sum = 0.0;
			for(int j=0;j<objNumber_;j++) {
				sum += sets[1].get(i).getNormalizedObjective(j)*lamada[j];
			}
			//sets[1].get(i).setFitness(sum);
			sets[1].get(i).setFitness(sets[1].get(i).getSumValue());
		}
		
	}
	
	/**
	 * Returns a <code>SolutionSet</code> containing the solutions of a given
	 * rank.
	 * 
	 * @param rank
	 *            The rank
	 * @return Object representing the <code>SolutionSet</code>.
	 */
	public SolutionSet getSubfront(int rank) {
		return ranking_.get(rank);
	} // getSubFront

	/**
	 * Returns the total number of subFronts founds.
	 */
	public int getNumberOfSubfronts() {
		return ranking_.size();
	} // getNumberOfSubfronts

}
