/**
 * dMOPSO_main.java
 * 
 * @version 1.0
 */

package jmetal.metaheuristics.dmopso;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.mutation.Mutation;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.UF1;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

public class dMOPSO_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args
	 *            Command line arguments. The first (optional) argument
	 *            specifies the problem to solve.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options -
	 *             jmetal.metaheuristics.mocell.MOCell_main -
	 *             jmetal.metaheuristics.mocell.MOCell_main problemName -
	 *             jmetal.metaheuristics.mocell.MOCell_main problemName
	 *             ParetoFrontFile
	 */
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
	      BufferedWriter bw      = new BufferedWriter(osw)        ;               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");
	        bw.newLine();        
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	public static void main(String[] args) throws JMException, IOException,
			ClassNotFoundException {
		

		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("dMOPSO_main.log");
		logger_.addHandler(fileHandler_);
		for(int fun=6;fun<=12;fun++){
			int runtimes=10;
			//double[] GDarray=new double[runtimes];
			  double[] IGDarray=new double[runtimes];
			  //double[] spreadarray=new double[runtimes];
			  //double[] Hypervolume=new double[runtimes];
			  long Execution_time=0;

		for(int i=0;i<runtimes;i++){

			Problem problem; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Mutation mutation; // "Turbulence" operator
			HashMap parameters; // Operator parameters
			QualityIndicator indicators; // Object to get quality indicators
		indicators = null;
		problem=null;
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
			      problem = new DTLZ1("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ1_M10.dat" ) ;
			    	}
		  	if(fun==7){
			      //problem = new DTLZ2("Real",10,2);
		  		problem = new DTLZ2("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==8){
			      problem = new DTLZ3("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==9){
			      problem = new DTLZ4("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}
		  	if(fun==10){
			      problem = new DTLZ5("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ5_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==11){
			      problem = new DTLZ6("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ5_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==12){
			      problem = new DTLZ7("Real",10,10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ7_M10.dat" ) ;
			    	}
			if(fun==13){
		  	      problem = new WFG1("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG1_605.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==14){
		  	      problem = new WFG2("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG2_111.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==15){
		  	      problem = new WFG3("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG3_right.txt" ) ;
		  	    	}
			if(fun==16){
			      problem = new WFG4("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG4_1181.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==17){
			      problem = new WFG5("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG5_694.txt" ) ;
			    	}
		  	if(fun==18){
		  	      problem = new WFG6("Real",4,8,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\wfg6_166.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==19){
			      problem = new WFG7("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG7_2435.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==20){
			      problem = new WFG8("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG8_201.txt" );
			    	}
			if(fun==21){
			      problem = new WFG9("Real",4,8,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\WFG9_2591.txt" ) ;
			    	}//problem = new WFG1("Real");
		  	if(fun==22){
			      problem = new Fonseca("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\Fonseca.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==23){
			      problem = new Kursawe("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\Kursawe.pf" ) ;
			    	}
			if(fun==24){
			      problem = new Schaffer("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\Schaffer.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==25){
			      problem = new UF1("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF1_500.txt" ) ;//.txt
			    	}
			if(fun==26){
			      problem = new UF2("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF2_500.txt" ) ;
			    	}
			if(fun==27){
			      problem = new UF3("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF3_500.txt" ) ;
			    	}
			if(fun==28){
			      problem = new UF4("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF4_500.txt" ) ;
			    	}
			if(fun==29){
			      problem = new UF5("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF5_21.txt" ) ;
			    	}
			if(fun==30){
			      problem = new UF6("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF6_668.txt" ) ;
			    	}
			if(fun==31){
			      problem = new UF7("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF7_500.txt" ) ;
			    	}
			if(fun==32){
			      problem = new UF8("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF8.DAT" ) ;
			    	}
			if(fun==33){
			      problem = new UF9("Real");
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF9.DAT" ) ;
			    	}
		} // else

		algorithm = new dMOPSO(problem);

		// Algorithm parameters
		algorithm.setInputParameter("swarmSize", 275);
		algorithm.setInputParameter("maxAge", 2);
		algorithm.setInputParameter("maxIterations", 100);
		algorithm.setInputParameter("functionType", "_PBI");
		
		algorithm.setInputParameter("dataDirectory",
				"E:\\sbMaOP\\SbMaOP\\Weights");

		parameters = new HashMap();

		// Execute the Algorithm
		long initTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		Execution_time += (System.currentTimeMillis() - initTime);

		// Result messages
		/*logger_.info("Total execution time: " + estimatedTime + "ms");
		logger_.info("Objectives values have been writen to file FUN");
		population.printObjectivesToFile("FUN");
		logger_.info("Variables values have been writen to file VAR");
		population.printVariablesToFile("VAR");

		if (indicators != null) {
			logger_.info("Quality indicators");
			logger_.info("Hypervolume: "
					+ indicators.getHypervolume(population));
			logger_.info("GD         : " + indicators.getGD(population));
			logger_.info("IGD        : " + indicators.getIGD(population));
			logger_.info("Spread     : " + indicators.getSpread(population));
			logger_.info("Epsilon    : " + indicators.getEpsilon(population));
		} // if*/
		//population.printObjectivesToFile("Run"+ i + "-FUN-" + problem.getName()
		//		+ "-DMOPSO-FUN");
		//GDarray[i]=indicators.getCECIGD(population);
	    IGDarray[i]=indicators.getIGD1(population);
	   // spreadarray[i]=indicators.getEpsilon(population);
	   // Hypervolume[i]=indicators.getHypervolume(population);
		}
	//printGD("FUN"+fun+"GD",GDarray);
		//printGD("FUN"+fun+"IGD",IGDarray);
		//printGD("FUN"+fun+"IGD",IGDarray);
		//printDiversity("FUN"+fun+"Diversity",spreadarray);
		//printHypervolume("FUN"+fun+"Hypervolume",Hypervolume);
		//double sumGD=0;
		double sumIGD=0;
		//double sumSP=0;
		//double sumHypervolume=0;
		for(int i=0;i<runtimes;i++){
		  //sumGD+=GDarray[i];
		  sumIGD+=IGDarray[i];
		  //sumSP+=spreadarray[i];
		 // sumHypervolume+=Hypervolume[i];
		}
		logger_.info("Total execution time: " + Execution_time + "ms");
		//System.out.println("avrGD-fun"+fun+" = "+sumGD/runtimes);
		//System.out.println("avrSP-fun"+fun+" = "+sumSP/runtimes);
		System.out.println("avrIGD-fun"+fun+"= "+sumIGD/runtimes);
		//System.out.println("avrHV-fun"+fun+" = "+sumHypervolume/runtimes);
		}
	} // main
} // dMOPSO_main