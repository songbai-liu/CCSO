package jmetal.metaheuristics.ccso;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;

public class EnvironmentalSelection {
	SolutionSet union_;
	SolutionSet leaders_;
	int swarmSize_;
	int objNumber_;
	
	public EnvironmentalSelection(SolutionSet union, SolutionSet leaders, int swarmSize, int objNumb) {
		this.union_ = union;
		this.leaders_ = leaders;
		this.swarmSize_ = swarmSize;
		this.objNumber_ = objNumb;
	}
	
	public SolutionSet eliteSelection() {
		SolutionSet[] subPopulation = null;
		subPopulation = new MostSimilarBasedSampling(union_, objNumber_).getIdealPointOrientedPopulation(swarmSize_);
		SolutionSet[] ReferenceSet = new LeastSimilarBasedSampling(subPopulation[0],objNumber_).getIdealPointOrientedPopulation(objNumber_);
    	subPopulation[0].clear();
    	SolutionSet subSets[] = new SolutionSet[swarmSize_];
    	for(int i=0;i<objNumber_;i++){
    		subSets[i] = new SolutionSet();
    		subSets[i].add(ReferenceSet[0].get(i));
    		subPopulation[0].add(ReferenceSet[0].get(i));
    	}
    	
		for(int i=0; i<swarmSize_-objNumber_;i++){
			  subSets[i+objNumber_] = new SolutionSet();
			  subSets[i+objNumber_].add(ReferenceSet[1].get(i));
			  subPopulation[0].add(ReferenceSet[1].get(i));
		}
		
		for(int i=0; i<subPopulation[1].size();i++){
			Solution s1 = subPopulation[1].get(i);
			double minDis = computeDistance(s1,subPopulation[0].get(0));
			int minIndex = 0;
			for(int j=1; j<subPopulation[0].size();j++){
				Solution s2 = subPopulation[0].get(j);
				double angle = computeDistance(s1,s2);
				if(angle<minDis){
				   minDis = angle;
				   minIndex = j;
				}
			}
			subSets[minIndex].add(s1);
		}
		
		for(int i=0; i<objNumber_;i++){
			double rd = PseudoRandom.randDouble();
			if(rd > 0.9){
				leaders_.add(new Solution(subSets[i].get(0)));
			}else{
				double minValue; 
				minValue = subSets[i].get(0).getSumValue();
				int minIndex = 0;
				for(int j=1; j<subSets[i].size();j++){
				   double value; 
				   value = subSets[i].get(j).getSumValue();
				   if(value < minValue){
					  minValue = value;
					  minIndex = j;
				   }
				}
				leaders_.add(new Solution(subSets[i].get(minIndex)));
			}
		}
		for(int i=objNumber_; i<swarmSize_;i++){
		    double minValue; 
		    minValue = subSets[i].get(0).getSumValue();
			int minIndex = 0;
			for(int j=1; j<subSets[i].size();j++){
			   double value; 
			   value = subSets[i].get(j).getSumValue();
			   if(value < minValue){
				  minValue = value;
				  minIndex = j;
			   }
			}
			leaders_.add(new Solution(subSets[i].get(minIndex)));
		}
    	
		return this.leaders_;
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
