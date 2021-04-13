package jmetal.metaheuristics.ccso;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;

public class LeastSimilarBasedSampling {
	SolutionSet solutionSet_;
	int objNumber_;
	public LeastSimilarBasedSampling(SolutionSet solutionSet, int number){
		this.solutionSet_ = solutionSet;
		this.objNumber_ = number;
	}
	
	public SolutionSet[] getIdealPointOrientedPopulation(int individualSize){
		int solutionSize = solutionSet_.size();
		
    /*subSets[0] preserves the reference vector solutions 
    * and the remain solutions are saved in subSets[1]*/
		SolutionSet subSet = new SolutionSet();
		
		for(int p=0;p<solutionSize;p++){
			solutionSet_.get(p).setRemove(false);
		}
     
		/*get one initial solution vectors, and saved them in subSets[0]*/
		int id_max1 = 0;
		int id_max2 = 1;
		double maxAng = computeDistance(solutionSet_.get(0), solutionSet_.get(1));
		for(int k=0; k<solutionSet_.size(); k++){
			for(int s=k+1;s<solutionSet_.size();s++){
				double ang = computeDistance(solutionSet_.get(k), solutionSet_.get(s));
				if(ang > maxAng){
					maxAng = ang;
					id_max1 = k;
					id_max2 = s;
				}
			}
		}//for k
		solutionSet_.get(id_max1).setRemove(true);
		solutionSet_.get(id_max2).setRemove(true);
		subSet.add(solutionSet_.get(id_max1));
		subSet.add(solutionSet_.get(id_max2));
		
		double[] angles = new double[solutionSet_.size()];
		int index[] = new int[solutionSet_.size()];
	 /*compute the angle between each not removed solution in solutionSet_ and subSet[0]*/
		for(int i=0;i<solutionSet_.size();i++){
			angles[i] = -1.0;
			Solution sol2 = solutionSet_.get(i);
			if(!sol2.isRemove()){
				double minAngle = computeDistance(sol2, subSet.get(0));
				int minIndex = 0;
				for(int j=1;j<subSet.size();j++){
					Solution sol3 = subSet.get(j);
					double angle = computeDistance(sol2, sol3);
					if(angle < minAngle){
						minAngle = angle;
						minIndex = j;
					}
				}//for j
				angles[i] = minAngle;
				index[i] = minIndex;
			}//if
		}//for i
		
		//int remain = individualSize - objNumber_;
		int remain = individualSize - 2;
		while(remain > 0){
		 /*find the current solution with the maximum angle to subSets[0]*/
			double maxAngle = -1.0e+30;
			int maxAngleID = -1;
			for(int a=0;a<solutionSet_.size();a++){
				if(!solutionSet_.get(a).isRemove()){
					if(angles[a] > maxAngle){
						maxAngle = angles[a];
						maxAngleID = a;
					}
				}
			}//for a
		/*maximum angle based addition*/
			solutionSet_.get(maxAngleID).setRemove(true);
			subSet.add(solutionSet_.get(maxAngleID));
		/*update angles*/
		    for(int b=0;b<solutionSet_.size();b++){
		    	Solution sol4 = solutionSet_.get(b);
		    	if(!sol4.isRemove()){
		    		Solution sol5 = solutionSet_.get(maxAngleID);
		    		double angle = computeDistance(sol4,sol5);
		    		if(angle < angles[b]){
		    			angles[b] = angle;
		    			index[b] = subSet.size()-1;
		    		}
		    	}
		    }
		    remain--;
		}//while
		
		SolutionSet[] subSets = new SolutionSet[2];
		subSets[0] = new SolutionSet();
		subSets[1] = new SolutionSet();
		for(int i=0;i<individualSize;i++){
			subSets[0].add(subSet.get(i));
		}
		for(int i=0;i<solutionSet_.size();i++){
			if(!solutionSet_.get(i).isRemove()){
				subSets[1].add(solutionSet_.get(i));
			}
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
