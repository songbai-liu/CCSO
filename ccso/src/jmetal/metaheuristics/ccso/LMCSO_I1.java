// MCCSO.java
//
// Author:
//     Songbai Liu <Songbai209@qq.com>
// Copyright (c) 2019 Songbai Liu
//
// This Program is free software: you can redistribute it and/or modify 
// it under the terms of the GNU Lesser General Public License as published
// by the free software foundation, either version 4 of license, or any 
// later version (at your option).

package jmetal.metaheuristics.ccso;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.logging.Level;
import java.util.logging.Logger;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.qualityIndicator.Hypervolume;
import jmetal.qualityIndicator.util.MetricsUtil;
import jmetal.util.Distance;
import jmetal.util.DistanceShifted;
import jmetal.util.JMException;
import jmetal.util.Permutation;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.CrowdingDistanceComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.comparators.RankComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.ranking.ThetaRanking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class LMCSO_I1 extends Algorithm {
	/**
	 * Stores the number of particles_ used
	 */
	private int swarmSize_;
	/**
	 * Stores the maximum number of iteration_
	 */
	private int maxIterations_;
	/**
	 * Stores the current number of iteration_
	 */
	private int iteration_;
	/**
	 * Stores the particles
	 */
	private SolutionSet particles_;
	private SolutionSet leaders_;
	private SolutionSet union_;
	
	private Comparator dominance_;
	private Distance distance = new Distance();
	private SolutionSet groupW_;
	private SolutionSet groupL_;
	
	private double[][] speed_;
	double theta_;
	double[] zideal_;
	double[][] lamada_;

	private Operator polynomialMutation_;
	
	private double deltaMax_[];
	private double deltaMin_[];
	
	private MetricsUtil utils_;
	private Hypervolume hv_;
	
	public LMCSO_I1(Problem problem) {
		super(problem);
		this.utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
		this.hv_ = new Hypervolume();
	}

	/**
	 * Runs of the MCCSO algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initParams();
		// -> Step 1 Create the initial particle swarm and evaluate
		for (int i = 0; i < swarmSize_; i++) {
			Solution particle = new Solution(problem_);
			problem_.evaluate(particle);
			problem_.evaluateConstraints(particle);
			particles_.add(particle);
		}
		for (int i = 0; i < swarmSize_; i++) {
			Solution particle = particles_.get(i);
			for(int n=0;n<problem_.getNumberOfVariables();n++) {
				speed_[i][n] = 0.0;
			}
		}	
		// -> Step4 Iterations
		while (iteration_ < maxIterations_) {
			for (int i = 0; i < particles_.size(); i++) {
				Solution particle = new Solution(particles_.get(i));
				particle.setindex(i+swarmSize_);
				leaders_.add(particle);
				for(int n=0;n<problem_.getNumberOfVariables();n++) {
					speed_[i+swarmSize_][n] = speed_[i][n];
				}
			}
			//int remain = swarmSize_/2;
			//int index = 0;
			//NondominatedRanking ranking1 = new NondominatedRanking(particles_);
			//SolutionSet frontw = ranking1.getSubfront(index);
			estimateIdealPoint(particles_);
			translateObjective(particles_);
			sortWithHV(particles_,particles_);
			for(int k = 0; k < swarmSize_/2; k++) {
				groupW_.add(particles_.get(k));
			}
			for(int k=swarmSize_/2;k<swarmSize_;k++) {
				groupL_.add(particles_.get(k));
			}
		
			theta_ = 1.0 + (4.0*iteration_)/maxIterations_;
			//distance.crowdingDistanceAssignment(particles_,problem_.getNumberOfObjectives());
			// -> Step7 Compute the Speed_ for the particles
			try {
				computeSpeed();
			}catch (IOException ex) {
				Logger.getLogger(LMCSO_I.class.getName()).log(Level.SEVERE, null,
						ex);
			}
			// -> Step7 Compute the new position for the particles
			computeNewPositions();
			/**
			 * Apply a mutation operator to particles in the swarm
			 * 
			 * @throws JMException
			 */
			for (int i = 0; i < particles_.size(); i++) {
				polynomialMutation_.execute(particles_.get(i));
			}
			// Evaluate the new particles_ in new positions
			for (int i = 0; i < particles_.size(); i++) {
				Solution particle = particles_.get(i);
				problem_.evaluate(particle);
				problem_.evaluateConstraints(particle);
				iteration_++;
			}
			union_ = ((SolutionSet) particles_).union(leaders_);
			SolutionSet st = getStSolutionSet(union_,swarmSize_);
			leaders_.clear();
			particles_.clear();
			groupW_.clear();
			groupL_.clear();
			
		    sortWithHV(st,st);
			for (int k = 0; k < swarmSize_; k++) {
				particles_.add(st.get(k));
			}
			for(int i=0;i<swarmSize_;i++){
				Solution particle = particles_.get(i);
				int id = particle.getindex();
				for(int n=0;n<problem_.getNumberOfVariables();n++) {
					speed_[i][n] = speed_[id][n];
				}
				particle.setindex(i);
			}
			union_.clear();	
			//iteration_++;
		}
		Ranking rank_Final = new NondominatedRanking(particles_);
		return rank_Final.getSubfront(0);
	}
	
	/**
	 * Initialize all parameter of the algorithm
	 */
	public void initParams() {
		swarmSize_ = ((Integer) getInputParameter("swarmSize")).intValue();
		maxIterations_ = ((Integer) getInputParameter("maxIterations")).intValue();
		polynomialMutation_ = operators_.get("mutation");
		iteration_ = 0;
		particles_ = new SolutionSet(swarmSize_);
		leaders_ = new SolutionSet(swarmSize_);
		union_ = new SolutionSet(2*swarmSize_);
		groupW_ = new SolutionSet(swarmSize_/2);
		groupL_ = new SolutionSet(swarmSize_/2);
		theta_ = 1.0;
		speed_ = new double[2*swarmSize_][problem_.getNumberOfVariables()];
		zideal_ = new double[problem_.getNumberOfObjectives()];
		deltaMax_ = new double[problem_.getNumberOfVariables()];
		deltaMin_ = new double[problem_.getNumberOfVariables()];
		for (int i = 0; i < problem_.getNumberOfVariables(); i++) {
			deltaMax_[i] = (problem_.getUpperLimit(i) - problem_.getLowerLimit(i)) / 2.0;
			deltaMin_[i] = -deltaMax_[i];
		} // for
		dominance_ = new DominanceComparator();
	}
	/*
	 * Estimate the Ideal Point 
	 */
	public void estimateIdealPoint(SolutionSet solutionSet){
		for(int i=0; i<problem_.getNumberOfObjectives();i++){
			zideal_[i] = 1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) < zideal_[i]){
					zideal_[i] = solutionSet.get(j).getObjective(i);
				}//if
			}//for
		}//for
	}
	
	public void translateObjective(SolutionSet solutionSet){
		double minNomal = Double.MAX_VALUE;
		int minIndex = 0;
		for(int i=0; i<solutionSet.size();i++){
			Solution sol = solutionSet.get(i);
			double sum = 0.0;
			double normal = 0.0;
			for(int j=0; j<problem_.getNumberOfObjectives();j++){
				double value = sol.getObjective(j) -  zideal_[j];
				sol.setNormalizedObjective(j, value);
				sol.setIthTranslatedObjective(j, value);
				sum += value;
				normal += value*value;
			}
			normal = Math.sqrt(normal);
			if(normal < minNomal) {
				minNomal = normal;
				minIndex = i;
			}
			sol.setDistanceToIdealPoint(normal);
			sol.setSumValue(sum);
		}
		double p = 1.0;
		for(int i=0; i<solutionSet.size();i++){
			double projection = 0.0;
			for(int j=0; j<problem_.getNumberOfObjectives();j++){
				projection += Math.pow(solutionSet.get(i).getIthTranslatedObjective(j), p);
			}
			projection = Math.pow(projection, 1/p);
			if(projection == 0) {
				projection = 0.00001;
			}
			for(int j=0; j<problem_.getNumberOfObjectives();j++){
				double value = solutionSet.get(i).getIthTranslatedObjective(j)/projection;
				solutionSet.get(i).setIthTranslatedObjective(j, value);
			}
		}
	}
	
	public void computeFitness(SolutionSet solutionSet) {
		DistanceShifted SDEdistance_  = new DistanceShifted();
		double [][] SDEdistance  =  SDEdistance_.translateDistanceMatrixShifted(solutionSet);
		int t = 1 ;
	  	for (int i = 0; i < SDEdistance.length; i++) {
	  		Arrays.sort(SDEdistance[i]);
	  		double tDistance =  SDEdistance[i][t] ; 
	  		solutionSet.get(i).setFitness(tDistance);                
	  	} // for
	}

	
	/**
	 * Update the speed of each particle
	 * 
	 * @throws JMException
	 */
	private void computeSpeed() throws JMException,IOException {
		int dim = problem_.getNumberOfVariables();
		int index1, index2;
		Solution winner1 = null;
		Solution winner2 = null;
		Solution loser = null;
		double fit1 = 0.0;
		double fit2 = 0.0;
		SolutionSet rivals = new SolutionSet();
		int size = groupW_.size();
		XReal teacher = null;
		XReal leader = null;
		for(int i=0;i<groupL_.size();i++) {
			loser = groupL_.get(i);
			loser.setLearningType(1);
			XReal particle = new XReal(loser);
			index1 =  loser.getindex();
			
			int id1 = PseudoRandom.randInt(0, groupW_.size()-1);
			int id2 = PseudoRandom.randInt(0, groupW_.size()-1); 
			while(id1 == id2) { 
				id2 = PseudoRandom.randInt(0, groupW_.size()-1); 
			}
			double rd = PseudoRandom.randDouble();
			if(rd <= 0.5) {
				fit1 = groupW_.get(id1).getDistanceToIdealPoint();
				fit2 = groupW_.get(id2).getDistanceToIdealPoint();
			}else {
				fit1 = -groupW_.get(id1).getCrowdingDistance();
				fit2 = -groupW_.get(id2).getCrowdingDistance();
			}
			if(fit1 < fit2) {
				winner1 = groupW_.get(id1);
				winner2 = groupW_.get(id2);
			}else {
				winner1 = groupW_.get(id2);
				winner2 = groupW_.get(id1);
			}
			teacher = new XReal(winner1);
			//leader = new XReal(winner2);
			
			 fit1 = getPBIDistance(loser,loser); 
			 for(int j=0;j<size;j++) { 
				 if(j != i){
					 Solution sol = groupW_.get(j); 
					 fit2 = getPBIDistance(sol, loser); 
					 if(fit2 < fit1) { 
						 rivals.add(sol); 
					 } 
				  } 
			 } 
			 int r = 0; 
			 if(rivals.size() > 0) { 
				 r = PseudoRandom.randInt(0, rivals.size()-1); 
				 winner1 = rivals.get(r); leader = new XReal(winner1); 
			}else { 
				r = PseudoRandom.randInt(0, size-1); 
				winner1 = groupW_.get(r); 
				leader = new XReal(winner1); 
			}
			double r1 = PseudoRandom.randDouble();
			double r2 = PseudoRandom.randDouble();
			double r3 = PseudoRandom.randDouble();
			double deltaSpeed;
			for(int var=0;var<dim;var++) {
				deltaSpeed = r1*speed_[index1][var];
				deltaSpeed += r2*(teacher.getValue(var)-particle.getValue(var)); 
				deltaSpeed += r3*(leader.getValue(var)-particle.getValue(var)); 
				if(deltaSpeed > deltaMax_[var]) { deltaSpeed = deltaMax_[var]; }
				if(deltaSpeed < deltaMin_[var]) { deltaSpeed = deltaMin_[var]; }
				speed_[index1][var] = deltaSpeed;
			}
			rivals.clear();
		}
		
		for(int i=0;i<size;i++) {
			groupW_.get(i).setLearningType(0);
		}
	}
	
	/**
	 * Update the position of each particle
	 * 
	 * @throws JMException
	 */
	private void computeNewPositions() throws JMException {
		int dim = problem_.getNumberOfVariables();
		for (int i = 0; i < swarmSize_; i++) {
			Solution sol = particles_.get(i);
			XReal particle = new XReal(sol);
			int learningType = sol.getLearningType();
			if(learningType != 0){
				for (int var = 0; var < dim; var++) {
					particle.setValue(var, particle.getValue(var) + speed_[i][var]);
					if (particle.getValue(var) < problem_.getLowerLimit(var)) {
						particle.setValue(var, problem_.getLowerLimit(var));
						speed_[i][var] = -1*speed_[i][var]; 
					}
					if (particle.getValue(var) > problem_.getUpperLimit(var)) {
						particle.setValue(var, problem_.getUpperLimit(var));
						speed_[i][var] = -1*speed_[i][var]; 
					}
				}
			}
		}
	}
	
	public double getPBIDistance(Solution sol, Solution ref) {
		double ip = 0;
		double refLenSQ = 0;
		double norm = 0.0;
		double distance = 0.0;
		double[] d = new double[2];
		for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
			ip += sol.getNormalizedObjective(j) * ref.getNormalizedObjective(j);
			refLenSQ += (ref.getNormalizedObjective(j) * ref.getNormalizedObjective(j));
			norm += sol.getNormalizedObjective(j) * sol.getNormalizedObjective(j);
		}
		refLenSQ = Math.sqrt(refLenSQ);
        norm = Math.sqrt(norm);
		d[0] = Math.abs(ip) / refLenSQ;

		d[1] = 0;
	    d[1] = norm*norm - d[0]*d[0];
		d[1] = Math.sqrt(d[1]);
		distance = d[0] + theta_*d[1];
		return distance;
	}
	
	public SolutionSet getStSolutionSet(SolutionSet ss,int size) {
		Ranking ranking = new NondominatedRanking(ss);
		int remain = size;
		int index = 0;
		SolutionSet front = null;
		SolutionSet mgPopulation = new SolutionSet();
		front = ranking.getSubfront(index);
		while ((remain > 0) && (remain >= front.size())) {

			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			} // for

			// Decrement remain
			remain = remain - front.size();

			// Obtain the next front
			index++;
			if (remain > 0) {
				front = ranking.getSubfront(index);
			} // if
		}
		if (remain > 0) { // front contains individuals to insert
			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			}
		}

		return mgPopulation;
	}
	
	/**
	   * Calculates how much hypervolume each point dominates exclusively. The points
	   * have to be transformed beforehand, to accommodate the assumptions of Zitzler's
	   * hypervolume code.
	   * @param front transformed objective values
	   * @return HV contributions
	   */
	  public double[] hvContributions(double[][] front) {
	    int numberOfObjectives = problem_.getNumberOfObjectives();
	    double[] contributions = new double[front.length];
	    double[][] frontSubset = new double[front.length - 1][front[0].length];
	    LinkedList<double[]> frontCopy = new LinkedList<double[]>();
	      Collections.addAll(frontCopy, front);
	    double[][] totalFront = frontCopy.toArray(frontSubset);
	    double totalVolume = hv_.calculateHypervolume(totalFront, totalFront.length, numberOfObjectives);
	    for (int i = 0; i < front.length; i++) {
	      double[] evaluatedPoint = frontCopy.remove(i);
	      frontSubset = frontCopy.toArray(frontSubset);
	      // STEP4. The hypervolume (control is passed to java version of Zitzler code)
	      double hv = hv_.calculateHypervolume(frontSubset, frontSubset.length, numberOfObjectives);
	      double contribution = totalVolume - hv;
	      contributions[i] = contribution;
	      // put point back
	      frontCopy.add(i, evaluatedPoint);
	    }
	    return contributions;
	  }
	 
	  public void sortWithHV(SolutionSet setFront, SolutionSet rSet) {
			for (int j = 0; j < setFront.size(); j++) {
				setFront.get(j).setCrowdingDistance(0.0);
		    }
			if (setFront.size() > 1) {
				 double[][] frontValues = setFront.writeObjectivesToMatrix();
			     int numberOfObjectives = problem_.getNumberOfObjectives();
			  // STEP 1. Obtain the maximum and minimum values of the Pareto front
			     double[] maximumValues = utils_.getMaximumValues(rSet.writeObjectivesToMatrix(), numberOfObjectives);
			     double[] minimumValues = utils_.getMinimumValues(rSet.writeObjectivesToMatrix(), numberOfObjectives);
			  // STEP 2. Get the normalized front
			     double[][] normalizedFront = utils_.getNormalizedFront(frontValues, maximumValues, minimumValues);
			  // compute offsets for reference point in normalized space
			     //double[] offsets = new double[maximumValues.length];
			     /*for (int i = 0; i < maximumValues.length; i++) {
			        offsets[i] = 10.0 / (maximumValues[i] - minimumValues[i]);
			     }*/
			  // STEP 3. Inverse the pareto front. This is needed because the original metric by Zitzler is for maximization problems
			     double[][] invertedFront = utils_.invertedFront(normalizedFront);
			  // shift away from origin, so that boundary points also get a contribution > 0
			     /*for (double[] point : invertedFront) {
			        for (int i = 0; i < point.length; i++) {
			          point[i] += offsets[i];
			        }
			     }*/
			  // STEP 4. calculate contributions and sort
			     double[] contributions = hvContributions(invertedFront);
			     for (int i = 0; i < contributions.length; i++) {
			       // contribution values are used analogously to crowding distance
			    	 setFront.get(i).setCrowdingDistance(contributions[i]);
			     }
			     setFront.sort(new CrowdingDistanceComparator());
			} 
	  }
	  

}
