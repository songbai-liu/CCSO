//  LMPSO_main.java
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

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.nsgaII.NSGADE;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.LSMOP.LSMOP1;
import jmetal.problems.LSMOP.LSMOP2;
import jmetal.problems.LSMOP.LSMOP3;
import jmetal.problems.LSMOP.LSMOP4;
import jmetal.problems.LSMOP.LSMOP5;
import jmetal.problems.LSMOP.LSMOP6;
import jmetal.problems.LSMOP.LSMOP7;
import jmetal.problems.LSMOP.LSMOP8;
import jmetal.problems.LSMOP.LSMOP9;
import jmetal.problems.M2M.MOP1;
import jmetal.problems.M2M.MOP2;
import jmetal.problems.M2M.MOP3;
import jmetal.problems.M2M.MOP4;
import jmetal.problems.M2M.MOP5;
import jmetal.problems.M2M.MOP6;
import jmetal.problems.M2M.MOP7;
import jmetal.problems.MaF.MaF1;
import jmetal.problems.MaF.MaF13;
import jmetal.problems.MaF.MaF2;
import jmetal.problems.MaF.MaF3;
import jmetal.problems.MaF.MaF4;
import jmetal.problems.MaF.MaF5_Convex;
import jmetal.problems.MaF.MaF6;
import jmetal.problems.MaF.MaF7;
import jmetal.problems.MaF.MaF8;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

/**
 * Class to configure and execute the NSGA-II algorithm.
 * 
 * Besides the classic NSGA-II, a steady-state version (ssNSGAII) is also
 * included (See: J.J. Durillo, A.J. Nebro, F. Luna and E. Alba "On the Effect
 * of the Steady-State Selection Scheme in Multi-Objective Genetic Algorithms"
 * 5th International Conference, EMO 2009, pp: 183-197. April 2009)
 */

public class CCSO_main_5obj {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args
	 *            Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName
	 *             paretoFrontFile
	 */
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");//写到缓冲区
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	public static void printave(String path,double aveIGD,double varianceIGD,double aveHypervolume,double varianceHV){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path) ;
	      OutputStreamWriter osw = new OutputStreamWriter(fos) ;
	      BufferedWriter bw      = new BufferedWriter(osw)  ;            
	            
	     // for (int i = 0; i < IGD.length; i++) {  
	        
	        bw.write(aveIGD+" ");
	        bw.newLine(); 
	        bw.write(varianceIGD+" ");
	        bw.newLine();
	        bw.write(aveHypervolume+" ");
	        bw.newLine();  
	        bw.write(varianceHV+" ");
	        bw.newLine();
	        /* Close the file */
		      bw.close();
		    }catch (IOException e) {
		      Configuration.logger_.severe("Error acceding to the file");
		      e.printStackTrace();
		    }       
		  } // printave
	
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAII_main.log");
		logger_.addHandler(fileHandler_);
		int n = 5000;
		int Gmax = n*15000;
		for(int fun=29;fun<=37;fun++){
			int runtimes=2;
			double[] GDarray=new double[runtimes];
			double[] IGDarray=new double[runtimes];
			double[] spreadarray=new double[runtimes];
			double[] Hypervolume=new double[runtimes];
			long Execution_time=0;
			Problem problem=null; // The problem to solve
		    Algorithm algorithm; // The algorithm to use
		    Operator crossover; // Crossover operator
		    Operator mutation; // Mutation operator
		    Operator selection; // Selection operator

		    HashMap parameters; // Operator parameters

		    QualityIndicator indicators; // Object to get quality indicators

		    indicators = null;
		    if (args.length == 1) {
		    	Object[] params = { "Real" };
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    } // if
		    else if (args.length == 2) {
		    	Object[] params = { "Real" };
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    	indicators = new QualityIndicator(problem, args[1]);
		    } // if
		    else { // Default problem
			if(fun==1){
		  	      problem = new ZDT1("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT1_501.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==2){
		  	      problem = new ZDT2("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT2_501.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==3){
		  	      problem = new ZDT3("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT3_269.txt" ) ;
		  	    	}
		  	if(fun==4){
			      problem = new ZDT4("Real",10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT4_501.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==5){
			      problem = new ZDT6("Real",10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT6_774.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==6){
				//problem = new DTLZ1("Real",10,2);
			      problem = new DTLZ1("Real",7,3);
			      
			      indicators = new QualityIndicator(problem);
			    	}
		  	if(fun==7){
			      //problem = new DTLZ2("Real",10,2);
		  		problem = new DTLZ2("Real",12,3);
			      
			      indicators = new QualityIndicator(problem);
			    	}//problem = new WFG1("Real");
			if(fun==8){
			      problem = new DTLZ3("Real",12,3);
			      
			      indicators = new QualityIndicator(problem);
			    	}//problem = new WFG1("Real");
			if(fun==9){
			      problem = new DTLZ4("Real",12,3);
			      
			      indicators = new QualityIndicator(problem);
			    	}
		  	if(fun==10){
			      problem = new DTLZ5("Real",19,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ5_10D.txt" ) ;
			      //indicators = new QualityIndicator(problem);
			    	}//problem = new WFG1("Real");
			if(fun==11){
			      problem = new DTLZ6("Real",19,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ6_10D.txt" ) ;
			      //indicators = new QualityIndicator(problem);
			    	}//problem = new WFG1("Real");
			if(fun==12){
			      problem = new DTLZ7("Real",27,8);
			      
			      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\DTLZ7.500000p.10D.pf.txt" ) ;
			      indicators = new QualityIndicator(problem);
			    	}
			if(fun==13){
		  	      problem = new WFG1("Real",4,20,3);
		  	      
		  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG2.500000p.5D.pf.txt" ) ;
		  	      indicators = new QualityIndicator(problem);
		  	    	}//problem = new WFG1("Real");
		  	if(fun==14){
		  	      problem = new WFG2("Real",4,20,3);
		  	      
		  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG2.500000p.5D.pf.txt" ) ;
		  	      indicators = new QualityIndicator(problem);
		  	    	}//problem = new WFG1("Real");
		  	if(fun==15){
		  	      problem = new WFG3("Real",4,20,3);
		  	      
		  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG3.500000p.5D.pf.txt" ) ;
		  	      indicators = new QualityIndicator(problem);
		  	    	}
			if(fun==16){
			      problem = new WFG4("Real",4,20,3);
			      
			      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
			      indicators = new QualityIndicator(problem);
			    	}//problem = new WFG1("Real");
			if(fun==17){
			      problem = new WFG5("Real",4,20,3);
			      
			      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
			      indicators = new QualityIndicator(problem);
			    	}
		  	if(fun==18){
		  	      problem = new WFG6("Real",4,20,3);
		  	      
		  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
		  	      indicators = new QualityIndicator(problem);
		  	    	}//problem = new WFG1("Real");
		  	if(fun==19){
			      problem = new WFG7("Real",4,20,3);
			      
			      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
			      indicators = new QualityIndicator(problem);
			    	}//problem = new WFG1("Real");
			if(fun==20){
			      problem = new WFG8("Real",4,20,3);
			      
			      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" );
			      indicators = new QualityIndicator(problem);
			    	}
			if(fun==21){
			      problem = new WFG9("Real",4,20,3);
			      
			      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
			      indicators = new QualityIndicator(problem);
			    	}//problem = new WFG1("Real");
		 	if(fun==22){
			      problem = new MaF1("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==23){
			      problem = new MaF2("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==24){
			      problem = new MaF3("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==25){
			      problem = new MaF4("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==26){
			      problem = new MaF5_Convex("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==27){
			      problem = new MaF6("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==28){
			      problem = new MaF7("Real",22,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==29){
			      problem = new LSMOP1("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP1_5D_20000.txt");
			    	}
			if(fun==30){
			      problem = new LSMOP2("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP2_5D_20000.txt");
			    	}
			if(fun==31){
			      problem = new LSMOP3("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP3_5D_20000.txt");
			    	}
			if(fun==32){
			      problem = new LSMOP4("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP4_5D_20000.txt");
			    	}
			if(fun==33){
			      problem = new LSMOP5("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP5_5D_20000.txt");
			    	}
			if(fun==34){
			      problem = new LSMOP6("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP6_5D_20000.txt");
			    	}
			if(fun==35){
			      problem = new LSMOP7("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP7_5D_20000.txt");
			    	}
			if(fun==36){
			      problem = new LSMOP8("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP8_5D_20000.txt");
			    	}
			if(fun==37){
			      problem = new LSMOP9("Real",n,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP9_5D_20000.txt");
			    	}
			if(fun==38){
			      problem = new UF1("Real",500);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF1_500.txt");
			    	}
			if(fun==39){
			      problem = new UF2("Real",500);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF2_500.txt");
			    	}
			if(fun==40){
			      problem = new UF3("Real",500);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF3_500.txt");
			    	}
			if(fun==41){
			      problem = new UF4("Real",500);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF4_500.txt");
			    	}
			if(fun==42){
			      problem = new UF5("Real",500,2,0.1);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF5_21.txt");
			    	}
			if(fun==43){
			      problem = new UF6("Real",500,2,0.1);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF6_668.txt");
			    	}
			if(fun==44){
			      problem = new UF7("Real",500);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF7_500.txt");
			    	}
			if(fun==45){
			      problem = new UF8("Real",100);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==46){
			      problem = new UF9("Real",100,0.1);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==47){
			      problem = new UF10("Real",100);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==48){
			      problem = new MOP1("Real");
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==49){
			      problem = new MOP2("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==50){
			      problem = new MOP3("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==51){
			      problem = new MOP4("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==52){
			      problem = new MOP5("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==53){
			      problem = new MOP6("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==54){
			      problem = new MOP7("Real");;
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			
		} // else
		for(int i=0;i<runtimes;i++){
		algorithm = new CCSOC(problem);
		//algorithm = new NSGADE(problem);
		// Algorithm parameters
		if(fun<=5){//ZDT1-ZDT6
			algorithm.setInputParameter("swarmSize", 100);
			algorithm.setInputParameter("maxIterations", 100);	
		} else if(fun<=12){//DTLZ1-DTLZ7
			algorithm.setInputParameter("swarmSize", 100);
			algorithm.setInputParameter("maxIterations", 1000);	
		}else if(fun<=21){//WFG1-WFG9
			algorithm.setInputParameter("swarmSize", 300);
			algorithm.setInputParameter("maxIterations", 500);
			algorithm.setInputParameter("T", 30);
		}else if(fun<=28){//MaF1-MaF7
			algorithm.setInputParameter("swarmSize", 2000);
			algorithm.setInputParameter("maxIterations", 2000);
		}else if(fun<=37){//LSMOP1-LSMOP9
			algorithm.setInputParameter("swarmSize", 336);
			algorithm.setInputParameter("maxIterations", Gmax);
			algorithm.setInputParameter("T", 20);
			algorithm.setInputParameter("div1", 5);
			algorithm.setInputParameter("div2", 6);
		}else if(fun<=44){//UF1-UF7
			algorithm.setInputParameter("swarmSize", 300);
			algorithm.setInputParameter("maxIterations", 5000);
			algorithm.setInputParameter("T", 20);
			algorithm.setInputParameter("div1", 299);
			algorithm.setInputParameter("div2", 0);
		}else if(fun<=47){//UF8-UF10
			algorithm.setInputParameter("swarmSize", 300);
			algorithm.setInputParameter("maxIterations", 5000);
			algorithm.setInputParameter("T", 30);
			algorithm.setInputParameter("div1", 23);
			algorithm.setInputParameter("div2", 0);
		}else if(fun<=52){//MOP1-MOP5
			algorithm.setInputParameter("swarmSize", 100);
			algorithm.setInputParameter("maxIterations", 3000);
			algorithm.setInputParameter("T", 20);
			algorithm.setInputParameter("div1", 99);
			algorithm.setInputParameter("div2", 0);
		}else{//MOP6-MOP7
			algorithm.setInputParameter("swarmSize", 300);
			algorithm.setInputParameter("maxIterations", 1000);
			algorithm.setInputParameter("T", 30);
			algorithm.setInputParameter("div1", 23);
			algorithm.setInputParameter("div2", 0);
		}
		
		
		
		// Mutation and Crossover for Real codification
		// Crossover operator
		parameters = new HashMap();
		parameters.put("CR", 1.0);
		parameters.put("F", 0.5);
		crossover = CrossoverFactory.getCrossoverOperator(
				"DifferentialEvolutionCrossover", parameters);		

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation",
				parameters);

		// Selection Operator
		parameters = null;
		selection = SelectionFactory.getSelectionOperator("BinaryTournament",
				parameters);

		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);

		// Add the indicator object to the algorithm
		algorithm.setInputParameter("indicators", indicators);

		// Execute the Algorithm
		long initTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		Execution_time+=(System.currentTimeMillis() - initTime);

		// Result messages
		population.printObjectivesToFile("CCSOC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1+8)+".txt" );
	    IGDarray[i]=indicators.getIGD1(population);
		}
	   /* wfghvCalculator1 wfg = new wfghvCalculator1(population,fun);
		  Hypervolume[i] = wfg.calculatewfghv();
		}
		printGD("NSGAII_10Obj_"+problem.getName()+"_HV.txt",Hypervolume);
		printGD("NSGAII_10Obj_"+problem.getName()+"_IGD.txt",IGDarray);*/

		  double sumHypervolume=0;
		  double sumIGD=0;
		  for(int i=0;i<runtimes;i++){
			  //sumHypervolume+=Hypervolume[i];
			  sumIGD+=IGDarray[i];
		  }	  	  
		  //aveHV   = sumHypervolume/runtimes; 
		 // printave("NSGAII_10Obj_Fun"+fun+".txt",aveIGD,varianceIGD,aveHV,varianceHV);
		  logger_.info("Total execution time: " + Execution_time + "ms");
		 // System.out.println("avrHV-fun"+fun+"= "+sumHypervolume/runtimes);
		  System.out.println("avrIGD-fun"+fun+"= "+sumIGD/runtimes);
	  } //main
	}
	} // NSGAII_main
