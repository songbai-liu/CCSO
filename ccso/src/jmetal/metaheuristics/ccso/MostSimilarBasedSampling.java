package jmetal.metaheuristics.ccso;

import java.util.ArrayList;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;
/*divid the population into two subSets based on toward the ideal point direction 
 * by the minimum-angle-delete method*/ 
public class MostSimilarBasedSampling {
	//List<SolutionSet> list = new <SolutionSet> ArrayList();
	SolutionSet solutionSet_;
	int objNumber_;
	
	public MostSimilarBasedSampling(SolutionSet solutionSet, int number){
		this.solutionSet_ = solutionSet;
		this.objNumber_ = number;
	}
	public SolutionSet[] getIdealPointOrientedPopulation(int individualSize){
		if(solutionSet_.size() < individualSize) {
			System.out.println("Error: size of the solutionSet in not enough");
		}
		int solutionSize = solutionSet_.size();
		for(int p=0;p<solutionSize;p++){
			solutionSet_.get(p).setRemove(false);
		}
		double minDist = Double.MAX_VALUE;
		double[] minDists = new double[solutionSize];//Angles between each solution and its closest solution 
		int[] minIndexs = new int[solutionSize];//Index of each solution's closest solution
		int[] index = new int[4];
		index[0] = index[1] = index[2] = index[3] = -1;
		
		for(int i=0;i<solutionSize;i++){
			boolean flag1 = false;
			double min = Double.MAX_VALUE;
			for(int j=0;j<solutionSize;j++){
				if(i != j){
					double dist = computeDistance(solutionSet_.get(i),solutionSet_.get(j));
					if(min > dist){
						flag1 = true;
						min = dist;
						index[0] = i;
						index[1] = j;
					}
				}
			}
			if(flag1 == false){
				System.out.println("Error: the distance between " + index[0] + " and " + index[1] + " is infinity");
			}
			minDists[i] = min;
			minIndexs[i] = index[1];
			if(minDists[i] < minDist){
				minDist = minDists[i];
				index[2] = index[0];
				index[3] = index[1];
			}
		}
		
		int[] secIndex = new int[2];
		double[] secMinDists = new double[solutionSize];
		int[] secMinIndexs = new int[solutionSize];
		secIndex[0] = -1;
		secIndex[1] = -1;
		for(int i=0;i<solutionSize;i++){
			boolean flagSec = false;
			double min = Double.MAX_VALUE;
			for(int j=0;j<solutionSize;j++){
				if(i != j && minIndexs[i] != j){
					double dist = computeDistance(solutionSet_.get(i),solutionSet_.get(j));
					if(min > dist){
						flagSec = true;
						min = dist;
						secIndex[0] = i;
						secIndex[1] = j;
					}
				}
			}
			if(flagSec == false){
				System.out.println("Error: the distance between " + secIndex[0] + " and " + secIndex[1] + " is infinity");
			}
			secMinDists[i] = min;
			secMinIndexs[i] = secIndex[1];
		}
		int remainSize = solutionSize;
		while(remainSize > individualSize){
			double fitness1;
			double fitness2;
			fitness1 = secMinDists[index[2]];
			fitness2 = secMinDists[index[3]];
			if(index[2] == index[3]){
				System.out.println("Error: The current two most similar solutions are the same solution " + index[2]);
			}
			if(fitness1 < fitness2){
				if(solutionSet_.get(index[2]).isRemove()){
				  System.out.println("Error: the solution " +index[2]+  " is deleted repeatedly");
				}
				solutionSet_.get(index[2]).setRemove(true);
				/*
				 * update the index of solution i whose minIndexs[i] = index[2]
				 */
				for(int i=0; i<solutionSize; i++){
					if(secMinIndexs[i] == index[2] && !solutionSet_.get(i).isRemove()){
						double sbmin = Double.MAX_VALUE;
						int sb = -1;
						//double sbAngle = 0.0;
						for(int j=0; j<solutionSize; j++){
							if(!solutionSet_.get(j).isRemove() && i!=j && minIndexs[i] != j){
								double sbDist = computeDistance(solutionSet_.get(i),solutionSet_.get(j));
								if(sbmin > sbDist){
									sbmin = sbDist;
									sb = j;
								}//if
							}
						}
						secMinDists[i] = sbmin;
						secMinIndexs[i] = sb;
					}
				}
				for(int i=0; i<solutionSize; i++){
					if(minIndexs[i] == index[2] && !solutionSet_.get(i).isRemove()){
						minDists[i] = secMinDists[i];
						minIndexs[i] = secMinIndexs[i];
						
						double sbmin = Double.MAX_VALUE;
						int sb = -1;
						//double sbAngle = 0.0;
						for(int j=0; j<solutionSize; j++){
							if(!solutionSet_.get(j).isRemove() && i!=j && minIndexs[i] != j){
								double sbDist = computeDistance(solutionSet_.get(i),solutionSet_.get(j));
								if(sbmin > sbDist){
									sbmin = sbDist;
									sb = j;
								}//if
							}
						}
						secMinDists[i] = sbmin;
						secMinIndexs[i] = sb;
					}
				}
			}else{
				
				if(solutionSet_.get(index[3]).isRemove()){
					 System.out.println("Error: the solution " +index[3]+  " is deleted repeatedly");
				}
				solutionSet_.get(index[3]).setRemove(true);
				/*
				 * update the index of solution i whose minIndexs[i] = index[3]
				 */
				for(int i=0; i<solutionSize; i++){
					if(secMinIndexs[i] == index[3] && !solutionSet_.get(i).isRemove()){
						double sbmin = Double.MAX_VALUE;
						int sb = -1;
				
						for(int j=0; j<solutionSize; j++){
							if(!solutionSet_.get(j).isRemove() && i!=j && minIndexs[i] != j){
								double sbDist = computeDistance(solutionSet_.get(i),solutionSet_.get(j));
								if(sbmin > sbDist){
									sbmin = sbDist;
									sb = j;
								}//if
							}
						}
						secMinDists[i] = sbmin;
						secMinIndexs[i] = sb;
					}
				}
				
				for(int i=0; i<solutionSize; i++){
					if(minIndexs[i] == index[3] && !solutionSet_.get(i).isRemove()){
						minDists[i] = secMinDists[i];
						minIndexs[i] = secMinIndexs[i];
						
						double sbmin = Double.MAX_VALUE;
						int sb = -1;
					
						for(int j=0; j<solutionSize; j++){
							if(!solutionSet_.get(j).isRemove() && i!=j && minIndexs[i] != j){
								double sbDist = computeDistance(solutionSet_.get(i),solutionSet_.get(j));
								if(sbmin > sbDist){
									sbmin = sbDist;
									sb = j;
								}//if
							}
						}
						secMinDists[i] = sbmin;
						secMinIndexs[i] = sb;
					}
				}
			}
			/*
			 * update the current two solution with the minimum similarity, i.e., update index[2] and index[3])
			 */
			double fDist = Double.MAX_VALUE;
			boolean flag3 = false;
			for(int p=0;p<solutionSize;p++){
				if(!solutionSet_.get(p).isRemove()){
					if(fDist > minDists[p]){
						flag3 = true;
						fDist = minDists[p];
						index[2] = p;
						index[3] = minIndexs[p];
					}
				}
			}
			if(flag3 == false){
				System.out.println("Error: the distance between " +index[3]+ " and " + index[3] + " is infinity");
			}
			if(solutionSet_.get(index[3]).isRemove()){
				 System.out.println("Error: the solution " +index[3]+  " will be deleted repeatedly");
			}
			
			remainSize = remainSize - 1;
			if(remainSize < individualSize) {
				System.out.println("Error: the remainSize is not enough");
			}
		}//while
		
		/*subSet1 preserves the reference vector solutions and the remain solutions are saved in subSet2*/
		SolutionSet[] subSets = new SolutionSet[2];
		subSets[0] = new SolutionSet();
		subSets[1] = new SolutionSet();
		for(int i=0;i<solutionSet_.size();i++){
			if(solutionSet_.get(i).isRemove()){//solutions not removed are reference vectors
				subSets[1].add(solutionSet_.get(i));
			}else{
				subSets[0].add(solutionSet_.get(i));
			}
		}
		if(subSets[0].size() != individualSize){
			System.out.println("Theoretical size = " + individualSize);
			System.out.println("Actual size = " + subSets[0].size());
			System.out.println("the other set size = " + subSets[1].size());
		}
		return subSets;
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
}
