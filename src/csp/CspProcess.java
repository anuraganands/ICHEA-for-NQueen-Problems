/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package csp;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.activation.UnsupportedDataTypeException;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import matlabcontrol.MatlabConnectionException;
import org.jdesktop.application.Application;
//import org.jfree.data.category.DefaultCategoryDataset;
//import org.jfree.data.xy.XYDataset;
//import org.jfree.data.xy.XYSeries;
//import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author Anurag
 */
public class CspProcess{
    //static boolean abToHoJaFlag;
    public boolean bMatlabDraw;
    private ArrayList<Chromosome> chromosomes_;
    private Queue<Chromosome> suspended_;
    public static Chromosome bestSoFar;
    private double prevBest_;
    private double curBest_;
    private int stillSameBestCount;
    private ArrayList<Chromosome> solutions_;
    private ArrayList<ArrayList<Double>> chromeValues;
    private UserInput userInput_;
    private ExternalData externalData_;
    //private int population_;
    //private int generation_;
    private int poolSize_;
    private int tourSize_;
    private int knearest_;
    private String dataType_;
    private Double[] range_;
    private final int MAX_MOVES = 10;
    private final int MUM = 20;
    private final int MU = 20;
    private double MUTATION_RATE = 0.1;
    //private final int ARCHIVE_MAX;
    private final double REPLACE_PERCENT = 0.10; //8% PERCENT of chromosomes replaced by new population
    private final double IMMUNITY_PERCENT = 0.10; 
    private final double PARTIAL_SOL_PERCENT = 0.10;
    private MyRandom r_;
    private ArrayList<ArrayList<Double>> sameBestChromeVals_; //stores top ranked SAME_BEST_VAL_PERCENT % of chromosomes
    private int hasAllSame_; //counter to check if the SAME_BEST_VAL_PERCENT % of chromosomes is same for SAME_BEST_GENERATIONS generations
    private final double SAME_BEST_VAL_PERCENT = 0.5; //top ranked percentage of total chrm population
    private final int SAME_BEST_GENERATIONS = 200;//50;//12//measure to SAME_BEST_VAL_PERCENT % of top ranked chromosomes remain same for number of generations.
    private final int NO_PROGRESS_LIMIT = 6;//limit for no progress made in the NO_PROGRESS_LIMIT generations.
    private boolean bStagnant;
    private int stagnantVisit;
    public static double bringCloserRatio = 0.5;
    private boolean bOptimizationMode;
    private double maxCSPval;
    public static ArrayList<ArrayList<ArrayList<Double>>> CSPsols;
    private final double FORCED_PERCENT = 0.75;
    public static boolean bInTransition;
    public static int dynamicConstraintNo;
    public static ArrayList<ArrayList<ArrayList>> dynamicConstraints = new ArrayList<ArrayList<ArrayList>> ();
    public static int MAX_FUNCTIONAL_CONSTRAINTS;
    private double tabuDist;
   static int MaxComb = 1000;//10; //
    private int MaxHospital = 1; //2; //MaxComb/2;
    public static int negFeasibleRange;
    private final int FIT_DP_LIMIT = 3; //3//6 //we use 3 for 10E-1 and 6 for 10E-3
    
    /**
     *
     * @param userInput
     * @throws MyException
     */
    public CspProcess(UserInput userInput) throws MyException{
        //this();
        this.userInput_ = userInput;
        this.externalData_ = null;

        if(userInput_ == null){
            throw new MyException("No user input provided.", "Incorrect Data",JOptionPane.ERROR_MESSAGE);
        }

        initialize();
    } // Toavoid calling this constructor

    public CspProcess(ExternalData externalData)throws MyException{
        //this();
        this.externalData_ = externalData;
        this.userInput_ = this.externalData_.getUserInput();

        if(userInput_ == null || this.externalData_ == null){
            throw new MyException("No user input provided or empty external data.", "Incorrect Data",JOptionPane.ERROR_MESSAGE);
        }

        initialize();
    }

    private CspProcess(){
        ;
    }

    /**
     * private Constructor used only for default initialization
     */
    private void initialize() throws MyException{
        //abToHoJaFlag = false;
        bMatlabDraw = false;
        chromosomes_ = new ArrayList<Chromosome>();
        suspended_ = new LinkedList<Chromosome>();
        solutions_ = new ArrayList<Chromosome>();
        this.tourSize_ = 2; //default value assumed
        this.knearest_ = (int)(0.05*userInput_.population); //default value assumed
        this.r_ = new MyRandom();
        
        this.poolSize_ = userInput_.population/2; //default values assumed
        //ARCHIVE_MAX = userInput_.population/2;
        this.dataType_ = this.userInput_.dataType;

        this.range_ = new Double[userInput_.totalDecisionVars];
        for (int i = 0; i < userInput_.totalDecisionVars; i++) {
            this.range_[i] = 0.5; //double assumed.
        }

        if (this.userInput_.population < 5 || this.userInput_.generation < 1){
            throw new MyException("poulation size should be > 5 and generation should be > 1", "Input Data Error!",JOptionPane.ERROR_MESSAGE);
        }
        hasAllSame_ = 0;
        sameBestChromeVals_ = null;
        bStagnant = false;
        prevBest_ = Double.POSITIVE_INFINITY;
        curBest_ = Double.POSITIVE_INFINITY;
        stillSameBestCount = 0;
        stagnantVisit = 0;
        bOptimizationMode = false;
        maxCSPval = Double.MAX_VALUE;
        MAX_FUNCTIONAL_CONSTRAINTS = userInput_.totalConstraints - userInput_.totalDecisionVars;
        //for dynammicConstrainNo: 0 means incremental and userInput_.totalConstraints - userInput_.totalDecisionVars
        //means concurrent
        dynamicConstraintNo = 0;//userInput_.totalConstraints - userInput_.totalDecisionVars; //0; 
        negFeasibleRange = 0;
        tabuDist = -1.0;
        
        CSPsols = new ArrayList<ArrayList<ArrayList<Double>>>();
        CSPsols.add(new ArrayList<ArrayList<Double>>());
        for (int i = 0; i < userInput_.totalConstraints; i++) {
            CSPsols.get(0).add(new ArrayList<Double>());
        }
        bInTransition = false;
    }
    
   

    /**
     * Starts the whole process
     */
    public void start(JProgressBar pb, boolean saveChromes, ByRef nextPrefSuggestion) throws MyException{
        ArrayList<Chromosome> parents;
        ArrayList<Chromosome> offspring;
        ArrayList<Chromosome> temp;
        ArrayList<Chromosome> CSPsolsDB = new ArrayList<Chromosome>();
        Chromosome tempChrom;
        double startTime = 0.0;
        double endTime = 0.0;
        int totalSaved;
        int g = 0; //generation
        Draw d = null;       
        PrintWriter runHistory = null;
        
        
        try{
            runHistory = new PrintWriter("broyden-e-1_run_history_run01.csv");
            if(bMatlabDraw){
                d = new Draw();                
                d.draw(matlabPlotBuildConstraints());  
            }
            
            
            
            
            
            
            //<< to measure feasible region space
//            dynamicConstraintNo = 1;
//            initializeChromosomes(this.chromosomes_, 10000, g);
//            int feasibleCount = 0;
//            int total = chromosomes_.size();            
//            for (Chromosome ch : this.chromosomes_) {
//                if(ch.getRankComponents().size() == 0){ //vios == 0
//                    feasibleCount++;
//                }
//            }                        
//            System.out.println("ro: "+feasibleCount*1.0/total);                    
//            Application.getInstance().exit();
            //>>
            
            
            
            
            
            initializeChromosomes(this.chromosomes_, userInput_.population, g);
            bestSoFar = this.chromosomes_.get(0);
                                    
            startTime = System.nanoTime();
            startTime = startTime/Math.pow(10, 9);

            for (g = 1; g <= userInput_.generation; g++) {
                    CSPsolsDB = new ArrayList<Chromosome>(); 
                    if(bestSoFar.isSolution()){
                        int cloneCount = 0;
                        long muRate;
   
                        if(externalData_ != null){
                            for (int i = 0; i < 1; i++) {
                                for (Chromosome ch : chromosomes_) {
                                    if(ch.isSolution()){ 
                                        tempChrom = (Chromosome)ch.clone();
                                        mutationSwap(tempChrom, 1,10); //0.1*userInput_.totalConstraints);//Math.ceil(2*Math.exp(-0.01*g)));
                                        CSPsolsDB.add(tempChrom);
                                    }
                                }
                                for (Chromosome ch : chromosomes_) {
                                    if(ch.isSolution()){ 
                                        tempChrom = (Chromosome)ch.clone();
                                        mutationGroupSwap(tempChrom);
                                        CSPsolsDB.add(tempChrom);
                                    }
                                }
                            }    
                        }
                    }
                
//                    
//                    // <editor-fold defaultstate="collapsed" desc="PMX for NQUEEN Code">
//                    //<<PMX one...
//                        final int popSize = userInput_.population;
//                        final int queenNum = userInput_.totalConstraints;
//                        double fitValProp[] = new double[popSize];
//                        
//        
//                        
//                        double fitTot = 0;
//                        for (int j = 0; j < popSize; j++) {
//                            fitTot += ((queenNum * (queenNum - 1)) / 2 ) - chromosomes_.get(j).getFitnessVal(0);                            
//                        }
//                        
//                        for (int j = 0; j < popSize; j++) {
//                            fitValProp[j] = (((queenNum * (queenNum - 1)) / 2 ) - chromosomes_.get(j).getFitnessVal(0))/fitTot;                            
//                        }
//                        
//
//                        //int[][] matingPool = new int [popSize][queenNum];
//                        ArrayList<Chromosome> matingPool = new ArrayList<Chromosome>();
//                        int i = 0;
//                        boolean done = false;
//                        while (done == false){
//                            double mother = Math.random() * popSize;
//                            int mothersPosition = 0;
//                            double father = Math.random() * popSize;
//                            int fathersPosition = 0;
//                            if(done == false){
//                                for(int j = 0; j < popSize; j++){
//                                    if ((mother > fitValProp[j])){
//                                        mothersPosition = j;
//                                        //for (int l = 0; l < queenNum; l++){
//                                            //matingPool[i][l] = queenPop[mothersPosition][l];
//                                            matingPool.add((Chromosome)chromosomes_.get(mothersPosition).clone());
//                                        //}
//                                        j = popSize; //use break yaar...
//                                        i++;
//                                    }
//                                }
//                            }
//                            if(i >= popSize)
//                                done = true;
//                            if(done == false){
//                                for (int k = 0; k < popSize; k++){
//                                    if ((father > fitValProp[k])){
//                                        fathersPosition = k;
//                                        if(fathersPosition == mothersPosition && mothersPosition != 0){
//                                            fathersPosition = 0;
//                                        }
//                                        if (fathersPosition == mothersPosition && mothersPosition == 0 ){
//                                            fathersPosition = 1;
//                                        }
//                //                        for (int m = 0; m < queenNum; m++){
//                //                            matingPool[i][m] = queenPop[fathersPosition][m];
//                //                        }
//                                        matingPool.add((Chromosome)chromosomes_.get(fathersPosition).clone());
//                                        k = popSize;
//                                        i++;
//                                    }
//                                }
//                            }
//                            if(i >= popSize)
//                                done = true;
//                        //System.out.println("Mom is: " + mothersPosition + "Dad is: " + fathersPosition);
//                        }
//                        chromosomes_.clear();
//                        parents.clear();
//                        for(i = 0; i < popSize; i++){
//                            //for(int j = 0; j < queenNum; j++){
//                                //queenPop[i][j] = matingPool[i][j];
//                           // }
//                            parents.add(matingPool.get(i));
//                        }   
//                        
//                        
//                    //>>PMX done....
//                    // </editor-fold>
//                    
                    ArrayList<Chromosome> elite;
                    final int eliteSize = 10; //(int)Math.ceil(userInput_.population*0.1); // new Double().intValue();//10% of population
                    elite = new ArrayList<Chromosome>();
                    
                    for (int i = 0; i < eliteSize; i++) {
                        elite.add((Chromosome)chromosomes_.get(i).clone());
                    }
                    
                    //<<option 1 
                    parents = new ArrayList<Chromosome>();
                    for (int i = 0; i < chromosomes_.size(); i++) { // pool size = 50% of population
                        parents.add((Chromosome)chromosomes_.get(i));
                    }
                    int i;
                    offspring = new ArrayList<Chromosome>();                   
                    chromosomes_.clear();                    
                    chromosomes_.addAll(elite);
                    
                    if(g<5000){
                        offspring.addAll(interRaceCrossover(parents));
                    }else{                    
                        for (int j = 0; j < parents.size()-1; j+=2) {
                            offspring.addAll(interRaceCrossoverInteger(new ArrayList<Chromosome>(parents.subList(j, j+2))));
                        }
                    }
                    //>>
                    
                    //<<option 2 
//                    parents = noveltyTournamentSelection(); //select best parents only. very slow progress...
//                    offspring = new ArrayList<Chromosome>();
//                    chromosomes_.clear();
//                    chromosomes_.addAll(elite);                    
//                    for (int j = 0; j < parents.size()-1; j+=2) {
//                        offspring.addAll(interRaceCrossoverInteger(new ArrayList<Chromosome>(parents.subList(j, j+2))));
//                    }
                    //>>
        
                    
                    //offspring = interRaceCrossover(parents); //crossover selected parents only  
                    
//                    System.out.println("@@@: "+offspring.size());
                    mutation(offspring);//mutating crossovered offspring only                                                                            
                    parents.clear();//no longer needed               
                                       
                    chromosomes_.addAll(offspring);//include all offspring into the chrm set.                        
                    chromosomes_.addAll(CSPsolsDB);
                    
                    
                    
                    
                    //<<Remove duplicate
//                    int szb = chromosomes_.size();                    
//                    ArrayList<String> lstVal = new ArrayList<String>();
//                    ArrayList<Chromosome> newChromosomes = new ArrayList<Chromosome>();                    
//                    for (Chromosome ch : chromosomes_) {
//                        if(!lstVal.contains(ch.getSatisfaction().toString())){
//                            lstVal.add(ch.getSatisfaction().toString());
//                            newChromosomes.add(ch);
//                        }                    
//                    }                    
//                    chromosomes_ = newChromosomes;                    
//                    int sza = chromosomes_.size();                    
//                    if(sza!= szb){
//                        sza = sza;
//                    }
                    //>> duplicate removed
                    
                    
                    
                    
                    temp = new ArrayList<Chromosome>();
                    initializeChromosomes(temp, userInput_.population-chromosomes_.size(), g);
                    chromosomes_.addAll(temp);
                    
                    sortAndReplace(g);//sort according to least violation first then on ro value (novelty)    
                    
                    System.out.println(g+",best: "+bestSoFar.getFitnessVal(0)+ ", curConsts: " + userInput_.total__updatedConstraints);
                    runHistory.println(g+","+bestSoFar.getFitnessVal(0)+"," + userInput_.total__updatedConstraints);
                    
//                    System.out.println("Gen: "+g); //+"best rank: " + bestSoFar.getRank() 
//                            + ", fitness: "+ bestSoFar.getFitnessVal(0)
//                            + ", rank list: " + bestSoFar.getSatisfaction() + ", vals: " + bestSoFar.getValsCopy());
//                    bestSoFar.tempSortBy = userInput_.solutionBy;
//                    System.out.println(bestSoFar);
                    
//                    System.out.println("top ones...");
//                    for (i = 0; i < userInput_.population; i++) {
//                        //System.out.print(new DecimalFormat  ("#.##").format(chromosomes_.get(i).getFitnessVal(0)) +", ");
//                        System.out.print(MyMath.roundN(chromosomes_.get(i).getFitnessVal(0),4)+", "); //+" - " + chromosomes_.get(i).getRankComponents().size()+", ");// +"["+ chromosomes_.get(i).vals_.get(0) + "], ");           
//                    }  

                
                if(bMatlabDraw){
                    d.draw(matlabPlotBuildGeneration(g));
                }
                
                if(pb != null)
                    pb.setValue(pb.getMinimum()+(pb.getMaximum()-pb.getMinimum())*g/(userInput_.generation));
            }            
        }catch(SolutionFoundException SFE){
            System.out.println("\nSolution found at generation " + (g));
            System.out.println("Reason: " + SFE.getMessage());
            runHistory.println("\n\nSolution found at generation " + (g));
            runHistory.println("Reason: " + SFE.getMessage());
            
            try {
                if(bMatlabDraw){
                    d.draw(matlabPlotBuildGeneration(g));
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
            g = userInput_.generation;
            pb.setValue(pb.getMinimum()+(pb.getMaximum()-pb.getMinimum())*g/(userInput_.generation));            
        }catch(MyException me){
            me.showMessageBox();
        }catch(UnsupportedDataTypeException udte){
            udte.printStackTrace();
        }catch (MatlabConnectionException mce) {
            mce.printStackTrace();           
        }
        catch(Exception e){
            e.printStackTrace();
            throw new MyException("Exception raised in Start Process", "Check Start Process()",JOptionPane.ERROR_MESSAGE);           
        }finally{
            endTime = System.nanoTime();
            endTime = endTime/Math.pow(10, 9);
            System.out.println("Process time(Sec): " + MyMath.roundN(endTime - startTime,2));
            System.out.println("total chromosomes: "+chromosomes_.size());
            int evals = 0;
            Chromosome.totalAppEvals = Chromosome.totalAppEvals/userInput_.totalConstraints;// equal sized chromes
            Chromosome.totalAppEvals = 2*Chromosome.totalAppEvals/(userInput_.totalConstraints+1); //each ICHEA eval is (N+1)/2 times faster than others, hence divided. ICHEA N -> Ohters N*(N+1)/2            
            Chromosome.totalRefEvals = Chromosome.totalRefEvals/userInput_.totalConstraints;// equal sized chromes            
            evals = Chromosome.totalAppEvals + Chromosome.totalRefEvals;
            
            System.out.println("Total Evaluations: "+ evals);
            
            runHistory.println("Process time(Sec): " + MyMath.roundN(endTime - startTime,2));
            runHistory.println("total chromosomes: "+chromosomes_.size());
            runHistory.println("Total Evaluations: "+ evals);
            
            System.err.flush();
            System.out.flush();
            //Thread.currentThread().sleep(100);//sleep for 1000 ms

            setSolution();
            
            if(this.solutions_.isEmpty()){                            
                System.out.println("No Solution Found :( ****************");                
                System.out.println("best chromosomes\n" + bestSoFar);
                runHistory.println("No Solution Found :( ****************"); 
                runHistory.println("\nbest chromosomes\n" + bestSoFar);
//                System.out.println("val size: "+bestSoFar.getValsCopy().size());
            }else{
                System.out.println("Solution Found");                
                System.out.println("best chromosomes\n" + bestSoFar);
                runHistory.println("\n\nSolution Found"); 
                runHistory.println("\nbest chromosomes\n" + bestSoFar);
            }  
            
            runHistory.close();
            if(externalData_ != null){ //external data is used.  
                nextPrefSuggestion.setVal(String.valueOf(externalData_.getNextPrefLimit()));
                String fileName = "partial_solutions_pref_"+externalData_.getCurPref()+".ichea";
                try {
                    FileOutputStream fos;
                    fos = new FileOutputStream(fileName);
                    ObjectOutputStream oos = new ObjectOutputStream(fos);

                    Set<Chromosome> s = new LinkedHashSet<Chromosome>(chromosomes_);
                    chromosomes_ = new ArrayList<Chromosome>(s); 
                    totalSaved = Math.min(chromosomes_.size(),(int)(PARTIAL_SOL_PERCENT*userInput_.population));
                    
                    chromosomes_ = new ArrayList<Chromosome>(chromosomes_.subList(0, totalSaved));
                    chromosomes_.add(0, bestSoFar);
                    for (Chromosome c : chromosomes_) {
                        try {
                            c.refreshFitness(maxCSPval);
                        } catch (SolutionFoundException ex) {
                            ;
                        }
                    }
                    
                    oos.writeObject(chromosomes_);//chromosomes_.subList(0, (int)(PARTIAL_SOL_PERCENT*userInput_.population)));
                    oos.flush();
                    oos.close();
                    System.out.println("["+ totalSaved + "] chromosomes of data successfully Saved to File ["+fileName+"].");
                } catch (FileNotFoundException fnfe) {
                    System.err.println("Serialize Error! File cannot be created.");
                } catch (IOException ioe){
                    ioe.printStackTrace();
                    System.err.println("Serialize Error! Cannot write to the file ["+fileName+"].");
                }                
            }
        }    
    }

    
     public ArrayList<Chromosome> getSolution(){
        return this.solutions_;
    }
    
    public ArrayList<Chromosome> getChromosomes(){
        return this.chromosomes_;
    }
    
    /**
     * It is the measure of novely. The larger the ro value the more the novelty
     * in the search space.
     * @param chrome
     * @return
     * @throws MyException
     * @throws UnsupportedDataTypeException 
     */
    double getRoValue(Chromosome chrome) throws MyException, UnsupportedDataTypeException{
        double ro;
        if(this.dataType_.contains("Integer")){
            ro = getIntegerRoValue(chrome);
        }
        else if(this.dataType_.contains("Double")){
            ro = getDoubleRoValue(chrome);
        }else{
            throw new UnsupportedDataTypeException("Only supports Integer and Double data type");
        }
        return ro;
    }
    
    /**
     * method determines ro value for nominal data types. the higher the better.
     * @param chrome
     * @return
     * @throws MyException 
     */
    double getIntegerRoValue(Chromosome chrome) throws MyException{
        double ro;
        ArrayList<Integer> validChromosomesIdx = new ArrayList<Integer>();
        Double []dist;
        int tempKnearest;
        
        if (chromosomes_.isEmpty()){
            throw new MyException("No chromosme population", "Variable Initialization Error",JOptionPane.ERROR_MESSAGE);
        }
 
        for (int i = 0; i < chromosomes_.size(); i++) {
            if (chromosomes_.get(i).getRank() != this.userInput_.total__updatedConstraints)
                validChromosomesIdx.add(i);            
        }

        dist = new Double[validChromosomesIdx.size()];
        for (int i = 0; i < validChromosomesIdx.size(); i++) {
            dist[i] = MyMath.norm(chrome.getValsCopy(), chromosomes_.get(validChromosomesIdx.get(i)).getValsCopy(), MyMath.DIST_DUPLICATE_MISMATCH);
//            NOTE: I am using "SQUARE of distance" instead of just distance
//            because I will be using variance for ro_min.
            dist[i] = Math.pow(dist[i], 2);

        }
        Arrays.sort(dist);        
        // x1 itself is included in this set which should have the value 0.
        if(dist.length<=knearest_){ //////@Danger code............................
            tempKnearest = dist.length-1;
        }else{
            tempKnearest = knearest_;
        }
        ro = (1.0/tempKnearest) * MyMath.sum(dist, 0, tempKnearest);//Note: should not be knearest-1 as x1 itself is also included

        ro = MyMath.roundN(ro, 0); //to reduce so much variations...
        return ro;                     
    }
    
    /**
     * ro value determines the rank of novelty. The higher value the better.
     * @param chrome
     * @return
     * @throws MyException 
     */
    double getDoubleRoValue(Chromosome chrome) throws MyException{
        double ro;
        int maxViolation = this.userInput_.totalConstraints;
        ArrayList<Integer> validChromosomesIdx = new ArrayList<Integer>();
        Double []dist;
        int tempKnearest;
        
//        if(chrome.getRank() == maxViolation-1){
//            return -1.0;
//        }
        
        if (chromosomes_.isEmpty()){
            throw new MyException("No chromosme population", "Variable Initialization Error",JOptionPane.ERROR_MESSAGE);
        }
 
        for (int i = 0; i < chromosomes_.size(); i++) {
            if (chromosomes_.get(i).getRank() != maxViolation)
                validChromosomesIdx.add(i);            
        }

        dist = new Double[validChromosomesIdx.size()];
        for (int i = 0; i < validChromosomesIdx.size(); i++) {
            dist[i] = MyMath.norm(chrome.getValsCopy(), chromosomes_.get(validChromosomesIdx.get(i)).getValsCopy(), MyMath.DIST_EUCLEADIAN);
//            NOTE: I am using "SQUARE of distance" instead of just distance
//            because I will be using variance for ro_min.
            dist[i] = Math.pow(dist[i], 2);

        }
        Arrays.sort(dist);        
        // x1 itself is included in this set which should have the value 0.
        if(dist.length<=knearest_){ //////@Danger code............................
            tempKnearest = dist.length-1;
        }else{
            tempKnearest = knearest_;
        }
        ro = Math.pow(1.0/tempKnearest, 2) * MyMath.sum(dist, 0, tempKnearest);//Note: should not be knearest-1 as x1 itself is also included

        ro = MyMath.roundN(ro, 0);
        return ro;
    }

    
    
    private ArrayList<String> matlabPlotBuildGeneration(int generation) throws Exception{
        ArrayList<String> MatlabCommands = new ArrayList<String>();
        ArrayList<Double> drawXdata;
        ArrayList<Double> drawYdata;
        
        if (userInput_.totalDecisionVars != 2){
            throw new Exception("Draw Error! Only 2D draw is allowed.");
        } 
        
                
        if(generation <= 1){
            MatlabCommands.add("hdata = plot(0,0,\'+b\');");
        }

        MatlabCommands.add("title(\'Gen: " + generation + "\');");                

        drawXdata = new ArrayList<Double>();
        drawYdata = new ArrayList<Double>();
        for (Chromosome c : chromosomes_) {
            drawXdata.add(c.getVals(0));//x axis;
            drawYdata.add(c.getVals(1));//y axis
        }
        MatlabCommands.add("set(hdata,\'XData\'," + drawXdata.toString() + ");"); 
        MatlabCommands.add("set(hdata,\'YData\'," + drawYdata.toString() + ");");
        MatlabCommands.add("refreshdata(hdata)");
        MatlabCommands.add("drawnow;");

        return MatlabCommands;
    }
    
    private ArrayList<String> matlabPlotBuildConstraints(){
        ArrayList<String> commands = new ArrayList<String>();
        
                commands.add("hold on;");
        commands.add("x = [-100:0.1:100];");
        commands.add("br = 10.0;");
        commands.add("sr = 9.9;");
        commands.add("y1p = sqrt(br^2 - x.^2);");
        commands.add("y1m = -sqrt(br^2 - x.^2);");
        commands.add("y2p = sqrt(sr^2 - x.^2);");
        commands.add("y2m = -sqrt(sr^2 - x.^2);");
        
        commands.add("y3p = sqrt(br^2 - (x+2*br - 0.1).^2);");
        commands.add("y3m = -sqrt(br^2 - (x+2*br - 0.1).^2);");
        commands.add("y4p = sqrt(sr^2 - (x+2*br - 0.1).^2);");
        commands.add("y4m = -sqrt(sr^2 - (x+2*br - 0.1).^2);");


        commands.add("fig1 = gcf;"); //get current figure or create figure fig1
        commands.add("axes1 = axes('Parent',fig1);"); //Create axes
        
        commands.add("ylim(axes1,[-100 100]);");
        commands.add("box(axes1,'on');");
        commands.add("hold(axes1,'all');");
        commands.add("plot(x,y1p);");
        commands.add("plot(x,y1m,'Parent',axes1);");
        commands.add("plot(x,y2p,'Parent',axes1);");
        commands.add("plot(x,y2m,'Parent',axes1);");
        commands.add("plot(x,y3p,'Parent',axes1);");
        commands.add("plot(x,y3m,'Parent',axes1);");
        commands.add("plot(x,y4p,'Parent',axes1);");
        commands.add("plot(x,y4m,'Parent',axes1);");
                
        
//        commands.add("hold on;");
//        commands.add("X1 = [-100:1:100];");
//        commands.add("YMatrix1 = -9*X1+6.0;");
//        commands.add("YMatrix2 = 9*X1 - 1.0;");
////        commands.add("XMatrix3 = 1;");
////        commands.add("XMatrix4 = 0;");
//        commands.add("YMatrix3 = -9*X1+7.0;");
//        commands.add("YMatrix4 = 9*X1 - 2.0;");
//
//        commands.add("fig1 = gcf;"); //get current figure or create figure fig1
//        commands.add("axes1 = axes('Parent',fig1);"); //Create axes
//        
//        commands.add("ylim(axes1,[-100 100]);");
//        commands.add("box(axes1,'on');");
//        commands.add("hold(axes1,'all');");
//        commands.add("plot(X1,YMatrix1);");
//        commands.add("plot(X1,YMatrix2,'Parent',axes1);");
////        commands.add("plot(XMatrix3,X1,'Parent',axes1);");
////        commands.add("plot(XMatrix4,X1,'Parent',axes1);");
//        commands.add("plot(X1,YMatrix3,'Parent',axes1);");
//        commands.add("plot(X1,YMatrix4,'Parent',axes1);");
        
        return commands;
    }
    
    /**
     * Remove stagnant best values from the population of chromosomes.
     * preprocess - chrm population should be sorted and contains unique
     * chromosomes.
     * @param PERCENT - what PERCENT of chromosomes to be checked.
     * @param sameBestVals - Arraylist of same best val
     * @param sameGens - same best val for how many generations
     * @return
     */
//    private boolean isStagnant(){    
//        boolean allsame;
//        boolean bstagnant;
//        bstagnant = false;
//        allsame = false;
//        ArrayList<ArrayList<Double>> diverse;
//
//        if(sameBestChromeVals_ != null)
//            diverse = (ArrayList<ArrayList<Double>>)sameBestChromeVals_.clone();
//        else{
//            diverse = new ArrayList<ArrayList<Double>>();
//            sameBestChromeVals_ = new ArrayList<ArrayList<Double>>();
//        }
//
//        for (int ofsp = 0; ofsp < (int)(SAME_BEST_VAL_PERCENT*userInput_.population); ofsp++) {
//            diverse.add(chromosomes_.get(ofsp).getValsCopy());
//        }
//        
//        HashSet<ArrayList<Double>> hashSet = new HashSet<ArrayList<Double>>(diverse);
//        diverse = new ArrayList<ArrayList<Double>>(hashSet);
//
//        if(diverse.size() == sameBestChromeVals_.size()){
//            allsame = true;
//        }else{            
//            sameBestChromeVals_ = new ArrayList<ArrayList<Double>>();
//            for (int ofsp = 0; ofsp < (int)(SAME_BEST_VAL_PERCENT*userInput_.population); ofsp++) {
//                sameBestChromeVals_.add(chromosomes_.get(ofsp).getValsCopy());
//            }            
//        }
//
//        if(allsame){
//            hasAllSame_++;
//        }else{
//            hasAllSame_ = 0;
//        }
//
//        if(hasAllSame_ >= SAME_BEST_GENERATIONS){
//            hasAllSame_ = 0;
//            sameBestChromeVals_ = null;
//            bstagnant = true;
//        }
//
//        return bstagnant;
//    }

    private void initializeChromosomes(ArrayList<Chromosome> chromosome, final int SIZE, final int gen) throws SolutionFoundException, Exception {
        boolean bInitialStage;
        if(SIZE<=0){
            chromosome = null;
            return;
        }
        
        if (externalData_ != null){
            if(gen<10){
                bInitialStage = true;
            }else{
                bInitialStage = false;
            }
                
            chromosome.addAll(externalData_.initializeExternalChrmosomes(SIZE));
        }else{
            initializeChromosomesRandomly(chromosome,SIZE);
        }
    }

    /**
     * Initializes Chromosomes with random values
     */
    private void initializeChromosomesRandomly(ArrayList<Chromosome> chromosome, final int SIZE) throws SolutionFoundException, Exception{
        Object rand = null;
        Chromosome tempChromosome;
               
        
        for (int i = 0; i < SIZE; i++) {
            tempChromosome = new Chromosome(this.userInput_.solutionBy, this.userInput_);
            for (int j = 0; j < userInput_.totalDecisionVars; j++) {
                if (userInput_.dataType.contains("Integer")){
                    rand = r_.randVal(userInput_.minVals.get(j).intValue(), userInput_.maxVals.get(j).intValue());
                }else if (userInput_.dataType.contains("Double")){
                    if(Math.random()<0.5){
                        rand = r_.randVal((Double)userInput_.minVals.get(j), (Double)userInput_.maxVals.get(j));
                    }else{
                        if(Math.random() < 0.5){
                            rand = userInput_.minVals.get(j);
                        }else{
                            rand = userInput_.maxVals.get(j);
                        }
                    }
                }
                else{
                    System.err.println("Incorrect use of data types");
                    System.exit(1);
                }
                tempChromosome.appendVal((Double)rand,maxCSPval);
            }
            chromosome.add(tempChromosome);
        }        
    }

//    private void setFitness(ArrayList<Chromosome> chrm, final int SIZE){
//        for (int ofsp = 0; ofsp < SIZE; ofsp++) {
//            //ObjectiveFunction.definition(chrm.get(ofsp));
//            chrm.get(ofsp).setObjectiveFunctionVars();
//        }
//    }

    /**
     * noveltyTournamentSelection() - Tournament selection based on novelty
     * of the chrm in the population.
     * @return Returns ArrayList<Chromosome> of parent selected population 
     */
    private ArrayList<Chromosome> noveltyTournamentSelection() throws MyException, UnsupportedDataTypeException{
        ArrayList<Chromosome> candidates = new ArrayList<Chromosome>();
        ArrayList<Chromosome> parents = new ArrayList<Chromosome>(); // shoud have this.pool sizse
        ArrayList<Integer> temp;
        double csize0, csize1;
        double ro0, ro1;
        int candidate0dominates;
        
        if(this.tourSize_ != 2){
            throw new MyException("Tour Should be 2", "Inappropriate Tour Size",JOptionPane.ERROR_MESSAGE);
        }
                 
        for (int p = 0; p < this.poolSize_; p++) {
            //select tourSize_ chromosomes k.e 2 chromosomes randomly from the population
            temp = MyRandom.randperm(0, chromosomes_.size()-1);
            candidates.clear();
            for (int t = 0; t < this.tourSize_; t++) {                
                candidates.add(chromosomes_.get(temp.get(t)));
            }
            temp = null;            

            try{
                    //<< you can move it bottom .....
                        ro0 = getRoValue(candidates.get(0)); //do not need to use getRo function, check sortnreplace function if it has already been set in tempRo property.
                        ro1 = getRoValue(candidates.get(1));

                        if (ro0 > ro1 || candidates.get(0).getValsCopy().size() == 1) // the larger the better
                            parents.add(candidates.get(0));
                        else if (ro1 > ro0 || candidates.get(1).getValsCopy().size() == 1)
                            parents.add(candidates.get(1));
                        else{
                        //>>..........................
                    Chromosome.tmpSortBy = Chromosome.BY_FITNESS;
                    csize0 = candidates.get(0).getRank();
                    csize1 = candidates.get(1).getRank();
                    Chromosome.tmpSortBy = userInput_.solutionBy;

                    if (csize0 < csize1){ // the lower the better
                        parents.add(candidates.get(0));                
                    }
                    else if (csize1 < csize0){
                        parents.add(candidates.get(1));                
                    }
                    else{  
                        //
                        //??here??    
                        //    
                        candidate0dominates = 0;

                        if(candidates.get(0).isStagnant(this.NO_PROGRESS_LIMIT)){
                            parents.add(candidates.get(1));
                        }else if(candidates.get(1).isStagnant(this.NO_PROGRESS_LIMIT)){
                            parents.add(candidates.get(0));
                        }else{
                            temp = MyRandom.randperm(0, 1);
                            parents.add(candidates.get(temp.get(0)));
                        }
                    }                                    
                }
                }catch(Exception e){
                e.printStackTrace();
            }
//            //>>
        }                
        
        return parents;
    }
    
    /**
     * IMPROPER METHOD... NEEDS CORRECTION.... Checks if the solution for CSP has been achieved
     * @return Returns ArrayList<Chromosome> of solution chromosomes.
     */
    private void setSolution(){
        //ArrayList<ArrayList<Double>> duplicates = new ArrayList<ArrayList<Double>>();
        chromeValues = new ArrayList<ArrayList<Double>>();
        int beforeSize, afterSize;
                
        for (Chromosome chromosome : this.chromosomes_) {
            if(chromosome.isSolution()){ // no violations
//                beforeSize = chromeValues.size();  
//                chromeValues.add(chrm.getValsCopy());
//
//                HashSet<ArrayList<Double>> hashSet = new HashSet<ArrayList<Double>>(chromeValues);
//                chromeValues = new ArrayList<ArrayList<Double>>(hashSet);            
//                afterSize = chromeValues.size();
//            
//                if(afterSize>beforeSize){
                    this.solutions_.add(chromosome);
//                }            
            }
        } 
        
        if(bestSoFar.isSolution()){
            this.solutions_.add(bestSoFar);
        }
        
    }
    
    public String printChromeValues(){
        String str;
        
        str = Integer.toString(chromeValues.size()) + "\n";
        str += Integer.toString(this.userInput_.totalConstraints) + "\n";
        for (int i = 0; i < chromeValues.size(); i++) {
            for (int j = 0; j < chromeValues.get(i).size(); j++) {
                str += chromeValues.get(i).get(j).toString() + " ";                
            }
            str += "\n";            
        }
        return str;
    }
    
//    private Chromosome notVals(Chromosome in){
//        Chromosome out;
//        
//        Double[] temp = new Double[userInput_.totalConstraints];
//        for (int i = 0; i < temp.length; i++) {
//            temp[i] = i*1.0;            
//        }
//        
//        for (Double d : in.getValsCopy()) {
//            temp[d.intValue()] = -1.0;
//        }
//        
//        ArrayList<Double> notVal = new ArrayList<Double>();
//        
//        
//        for (int i = 0; i < temp.length; i++) {
//            if(temp[i]!=-1.0){
//                notVal.add(temp[i]);
//            }       
//        }
//        
//        out = (Chromosome)in.clone();
//        out.setVals(notVal, maxCSPval);
//        return out;
//    }

     /**
     * inter race crossover - offers crossover between 2 different constraint regions only
     * the offspring will have better or same constraint violation than their parents.
     * This process requires 2 parents that produce 2 offspring
     * @param parents list of parents
     * @return returns offspring
     * @throws MyException
     * @throws UnsupportedDataTypeException
     */
    private ArrayList<Chromosome> interRaceCrossover(final ArrayList<Chromosome> parents) throws MyException, UnsupportedDataTypeException,  SolutionFoundException{
        ArrayList<Chromosome> candidates = new ArrayList<Chromosome>(this.tourSize_);
        ArrayList<Chromosome> offspring = new ArrayList<Chromosome>();
        ArrayList<Integer> tempIntAL; 

        //ArrayList<Double> directions;
        //ArrayList<Double> approachDist = new ArrayList<Double>(1);
        
        //double maxDist;
        //double ratio;
        int count;
        
        if(this.tourSize_ != 2){
            throw new MyException("Tour Should be 2", "Inappropriate Tour Size",JOptionPane.ERROR_MESSAGE);
        }
        
        if(parents.isEmpty()){
            System.out.println("Sigh! no parents!");
        }
        
        for (int i = 0; i < userInput_.population/2; i++) {
            if(Math.random() < 1.0){
                //Randomly pick two parents.
                tempIntAL = MyRandom.randperm(0, parents.size()-1);
                candidates.clear();
                for (int t = 0; t < this.tourSize_; t++) {                
                    candidates.add((Chromosome)parents.get(tempIntAL.get(t)));
                }
                
                try{
                    //Note here we can make integer and double combined problem
                    //set as well.
                    if(dataType_.contains("Integer")){
                        count = this.tourSize_;
                        //while both parents belong to same Constraint region
                        //Drawback - this will case very few crossover + very few
                        //final solutions. That might affect the optimization
                        //problem where we need many candidate solutions.                 

//                        while(candidates.get(0).hasSameRankComponent(candidates.get(1))){// &&
//                                //candidates.get(0).getRank() == candidates.get(1).getRank()){
//                            candidates.remove(1);
//                            candidates.add(parents.get(tempIntAL.get(count)));
//                            count++;
//                            if(count >= parents.size()){
//                                break;                                    
//                            }
//                        }                        
                        offspring.addAll(interRaceCrossoverInteger(candidates));//only 1 move
                        
                        
                        
                    }else if(dataType_.contains("Double")){
                        //further filter for boundary intersections... 
                        count = this.tourSize_;

                        while(!candidates.get(0).isMarriageCompatible(candidates.get(1))){
                            candidates.remove(1);
                            candidates.add((Chromosome)parents.get(tempIntAL.get(count)).clone());
                            count++;
                            if(count >= parents.size()){                            
                                break;                            
                            }
                        }
                        
                        
                        
                        offspring.addAll(interRaceCrossoverDouble(this.MAX_MOVES, candidates));
						//NOTE: if input var size is very big you may want to use the following instead of top line.
                        //if(userInput_.totalDecisionVars < 10)
                            //offspring.addAll(interRaceCrossoverDoubleStagnant(this.MAX_MOVES, candidates));
                        //else
                            //offspring.addAll(interRaceCrossoverDouble(this.MAX_MOVES, candidates));
                        
                        
                        
                        
                    }else{
                        throw new UnsupportedDataTypeException("Only supports Integer and Double");
                    }         

                }catch (UnsupportedDataTypeException udte) {
                    throw new UnsupportedDataTypeException("Check your data type");
                }
//                catch (MyException me){
//                    me.printMessage();
//                }
            }
        }
        
        return offspring;
    }

    /**
     * interRaceCrossoverInteger - is used only with nominal data types. for integer data
     * use interRaceCrossoverDouble. it virtually moves 2 parents.
     * @param move
     * @param candidates - parents from which offspring are sought.
     * @return returns offspring from given candidate parents
     */
    private ArrayList<Chromosome> interRaceCrossoverInteger(final ArrayList<Chromosome> candidates) throws SolutionFoundException{
        ArrayList<Chromosome> offspring = new ArrayList<Chromosome>();
        ArrayList<Integer> idx = new ArrayList<Integer>();
        ArrayList<Double> noGoods;
        Chromosome tempChrome;
        int move;
        boolean isSol;
        
        if (candidates.size() != 2){
            throw new UnsupportedOperationException("Require only 2 parents");
        }
        
        isSol = false;
        for (Chromosome can : candidates) {
            if(can.isSolution()){
                isSol = true;
                break;
            }
        }
        
        //Check common values in both candidate parents.
        int [] commonVals = new int[userInput_.totalConstraints]; //all initialized to 0
        int constVal = 1;
        
        
        for (int j = 0; j < candidates.size(); j++) {
            for (double v : candidates.get(j).getValsCopy()) {
                commonVals[(int)v] += constVal;
            }
            constVal = constVal*10;
        }
        //Those commonVals that have value 11 as element value, it means that is common in both parents.
        //otherwise it will have 1 or 10 respectively for both parents.
       
        int prevLength, newLength;
        //Technique 1 - Append chromosomes- multi-offpring (0-n) afrom 2 parents. - 
        //<< Build up structure for satisfaction list
//        if(Math.random() < 0.5){//5 && !isSol){ //!bStagnant){ //obviously Math.random is always [0 1)
//        if(!bOptimizationMode){
        final int maxSz = (int)Math.ceil(userInput_.population*0.9);
        if(Math.random()<.90){//!isSol){candidates.get(0).getValsCopy().size()<maxSz || candidates.get(1).getValsCopy().size()<maxSz){ //
            constVal = 1;
            for (int j = 0; j < candidates.size(); j++) {
                tempChrome = (Chromosome)candidates.get((j+1)%tourSize_);
                prevLength = tempChrome.getValsCopy().size();
                
                for (int i = 0; i < commonVals.length; i++) {
                    if(commonVals[i]==constVal){
                        tempChrome.appendVal(i, maxCSPval);//NOTE ofsp want getSatisfaction value but in this case both are same
                    }                    
                }

                newLength = tempChrome.getValsCopy().size();
                
//                if(newLength < prevLength){ //val added
//                    noGoods = tempChrome.findNoGoods();
//                    for (int i = 0; i < noGoods.size(); i++) {                    
//                        tempChrome.appendVal(noGoods.get(i), maxCSPval);
//                    }
//                }
                offspring.add((Chromosome)tempChrome);
                constVal = constVal*10;
            }
        }
        else{
            
////            ArrayList<Double> vals;
////            ArrayList<Integer> pos;
//////            int temp;
////            
////            for (int j = 0; j < candidates.size(); j++) {
////                tempChrome = (Chromosome)candidates.get((j+1)%tourSize_).clone();
////                noGoods = tempChrome.findNoGoods();
////                vals = tempChrome.getValsCopy();
////                
////                if(vals.size()>userInput_.population/2){  
////                    pos = MyRandom.randperm(1, vals.size());
////                    for (int i = 0; i < noGoods.size(); i++) {                    
////                        tempChrome.appendVal(noGoods.get(i), maxCSPval);
////                    }
////                }
////                
//////                if(vals.size()>userInput_.population/2){                
//////                    pos = MyRandom.randperm(0, vals.size());
//////                    for (int i = 0; i < noGoods.size(); i++) {                    
//////                        vals.add(pos.get(i)+i,noGoods.get(i));
//////                    }
//////                    tempChrome.setVals(vals, maxCSPval);
//////                    offspring.add(tempChrome);
//////                }
////            }  
                    
            
            
            Chromosome p0L, p0R, p1L, p1R;
            //Chromosome init_p0;
            int part0, part1;
            final double PERCENT = Math.random() ;//0.90;

            try{
                part0 = (int)Math.ceil(candidates.get(0).getValsCopy().size()*PERCENT);
                part1 = (int)Math.ceil(candidates.get(1).getValsCopy().size()*PERCENT);

                p0L = (Chromosome)candidates.get(0).clone();
                p0R = (Chromosome)candidates.get(0).clone();
                //init_p0 = (Chromosome)p0.clone();
                p1L = (Chromosome)candidates.get(1).clone();
                p1R = (Chromosome)candidates.get(1).clone();
                
                if(part0 <=1 || part1 <=1){
                    p0L.appendVal(candidates.get(1).getVals(0), maxCSPval);
                    offspring.add(p0L);
                    p1L.appendVal(candidates.get(0).getVals(0), maxCSPval);
                    offspring.add(p1L);
                    return offspring;
                }               
                
                p0L.restructure(1-PERCENT, true, maxCSPval);
                p0R.restructure(PERCENT, true, maxCSPval);
                //init_p0 = (Chromosome)p0.clone();
                p1L.restructure(1-PERCENT, false, maxCSPval);
                p1R.restructure(PERCENT, false, maxCSPval);
                
                if(p0L.getValsCopy().isEmpty() || p1R.getValsCopy().isEmpty()){
                    System.out.println("stope reee...");
                }
                
                //add half one at a time
                for (int i = 0; i < p1L.getValsCopy().size(); i++) {                    
                    p0L.appendVal(p1L.getVals(i),maxCSPval); //here val and satisfaction values are same so getVal() method is fine here                    
                }
                
                //add half one at a time
                for (int i = 0; i < p0R.getValsCopy().size(); i++) {                    
                    p1R.appendVal(p0R.getVals(i),maxCSPval); //here val and satisfaction values are same so getVal() method is fine here                    
                }
                
                //rank check??????????????
//                if(p0.getRank() <= candidates.get(0).getRank())
                    offspring.add(p0L);
//                if(p1.getRank() <= candidates.get(1).getRank())
                    offspring.add(p1R);

            }catch(IndexOutOfBoundsException e){
                System.err.println("No offspring\n" + e.getLocalizedMessage());
            }         
        }

        
        return offspring;
    }
   
    /**
     * highest valid index of CSPsols.
     * It is possible that the index = highest+1 is in buildup process and 
     * not available for usage.
     * 
     * @return highest valid index. -1 means no valid index available
     */
    private int maxIdxAvailCSPsols(){        
        int highestIndx = CSPsols.size()-1; //can be 1 initially.        
        
        for (ArrayList<Double> valGrp: CSPsols.get(highestIndx)) {
            if(valGrp.isEmpty()){
                highestIndx--; //it can become -1
                break;
            }
        }
        
        //last index can get cleared in CCSPfns so we take the second last
        
        if(highestIndx>0){
            highestIndx--;
        }
        
        return highestIndx;
    }
        
    
    /**
     * interRaceCrossoverDouble - can be used for interger or double data types.
     * double crossover reqires 2 parents and generate 2 offspring
     * process: the original genes of parents are moved closer to each other until
     * the better or same ofsp.e. (less or equal violations) is reached. the number
     * of moves is determined by  move parameter
     * @param move - number of maximum moves until the better/same solution is reached
     * @param candidates Two parents
     * @return Two offspring
     * @throws UnsupportedDataTypeException
     */
    
    private ArrayList<Chromosome> interRaceCrossoverDouble(final int move, final ArrayList<Chromosome> candidates) throws SolutionFoundException, UnsupportedDataTypeException{
        ArrayList<Double> delta;
        ArrayList<Double> newDelta;
        ArrayList<Chromosome> offspring = new ArrayList<Chromosome>();
        Chromosome childChrome = null;

        if (candidates.size() != 2){
            throw new UnsupportedOperationException("Require only 2 parents");
        }
      
        Chromosome p1 = null, p2=null; //parent 1 and parent 2;
        int pickedCons = -1;
        int randPickIdx = 0; //0 is always available... worst case can be empty..
        final int highestValidCSPidx = maxIdxAvailCSPsols();
        int k;
        int trials = 0;
        boolean bMoved = false;  
        ArrayList<Double> initDelta;
        ArrayList<Double> neighborDelta;
        ArrayList<Chromosome> hospital = new ArrayList<Chromosome>();

        ArrayList<String> permutes;  
        ArrayList<ArrayList<Double>> combinations;
        double dtemp;
        boolean isInvalid = false;
        
        
        for (int j = 0; j < this.tourSize_; j++) {            
            p1 = candidates.get(j);
            if(highestValidCSPidx >= 0 && Math.random()<FORCED_PERCENT && !p1.isSolution() && p1.getRankComponents().size()>0){                                               
                pickedCons = MyRandom.randperm(0, p1.getRankComponents().size()-1).get(0);
                pickedCons = p1.getRankComponents().get(pickedCons);
                
                p2 = new Chromosome(this.userInput_.solutionBy,  userInput_);
                randPickIdx = MyRandom.randperm((int)Math.floor(0.5*highestValidCSPidx), highestValidCSPidx).get(0);
                try{ //
                    p2.setVals(CSPsols.get(randPickIdx).get(pickedCons), maxCSPval);   
                }catch(IndexOutOfBoundsException iobe){ //temp measure... see below for reason
                    //happens when dynmic tabu constraints are introduced, and CSPsols repertoir does not carry the solutions for these new constraints    
                    p2 = null;
                }
                
            }else{
                if(highestValidCSPidx == -1){ //I use then when it is difficult to find CSP
                    for (int i = 0; i < p1.getRankComponents().size(); i++) {
                        pickedCons = p1.getRankComponents().get(i);                        
                        
                        if(CSPsols.get(0).get(pickedCons).isEmpty()){
                            p2 = null;
                            continue;
                        }else{
                            p2 = new Chromosome(this.userInput_.solutionBy,  userInput_);
                            p2.setVals(CSPsols.get(0).get(pickedCons), maxCSPval);
                            break;
                        }
                    }
                }                                               
                if(p2 == null)
                    p2 = candidates.get((j+1)%tourSize_);  
            }            
                       
            if(p2 == null)
                p2 = candidates.get((j+1)%tourSize_);  
            
            initDelta = new ArrayList<Double>(MyMath.vectorSubtraction(p2.getValsCopy(), p1.getValsCopy()));
            hospital.clear();
            combinations = MyMath.getXrandomBinPermutes(userInput_.totalDecisionVars, MaxComb);
            isInvalid = false;
            for(ArrayList<Double> dims: combinations){
                isInvalid = false;
                                
                k = 0;

                neighborDelta = MyMath.vectorMultiplication(true, dims, initDelta);
                
                p2 = new Chromosome(this.userInput_.solutionBy, this.userInput_);
                p2.setVals(MyMath.vectorAddition(p1.getValsCopy(), neighborDelta), maxCSPval);

                delta = neighborDelta; //MyMath.vectorSubtraction(p2.getValsCopy(), p1.getValsCopy());
                
                ArrayList<Double> prevVals = null, newVals=null;
                int forceFind = 0;
                do{
                    for (k = 1; k <= move; k++){ //) move; k++) {
                        //find which direction to move?
                        childChrome = new Chromosome(this.userInput_.solutionBy, this.userInput_);                                        
                        newDelta = MyMath.constMultiplicationToVector(Math.pow(bringCloserRatio,k), delta); 
                        newVals = MyMath.vectorAddition(p1.getValsCopy(), newDelta);
                        childChrome.setVals(newVals, maxCSPval);
                   
                        bMoved = false;

                        if(!childChrome.myParent(p1)){ //there is a gap/black hole so no point searching further.
                            hospital.add(childChrome);
                            bMoved = true;
                            break;
                        }
                        if(childChrome.isSolution() || childChrome.getRank()<=p1.getRank() || bStagnant){
                            hospital.add(childChrome);
                            if(childChrome.isSolution()){
                                bestSoFar = childChrome;//.clone()                                    
                                throw new SolutionFoundException("Sol found during crossover...");
                            }
                            bMoved = true;
                            break;
                        }
                        prevVals = newVals;
                    }
                    if(bMoved && k>=1){
                        p1  = childChrome;
                        delta = MyMath.vectorSubtraction(p2.getValsCopy(), p1.getValsCopy());
                    }else{
                        forceFind = 1;
                    }
                    forceFind++;
                }while(forceFind < 1);
            }
            
            if(!hospital.isEmpty()){
                Collections.sort(hospital);                
                offspring.add(hospital.get(0)); //(Chromosome)hospital.get(0).clone());
                hospital.clear();
            }            
        }

        
        return offspring;
    }
    
    
//    private ArrayList<Chromosome> interRaceCrossoverDouble(final int move, final ArrayList<Chromosome> candidates) throws SolutionFoundException, UnsupportedDataTypeException{
//        ArrayList<Double> delta;
//        ArrayList<Double> newDelta;
//        ArrayList<Chromosome> offspring = new ArrayList<Chromosome>(tourSize_);
//        Chromosome childChrome = null;
//
//        if (candidates.size() != 2){
//            throw new UnsupportedOperationException("Require only 2 parents");
//        }
//      
//        Chromosome p1 = null, p2=null; //parent 1 and parent 2;
//        int pickedCons = -1;
//        int randPickIdx = 0; //0 is always available... worst case can be empty..
//        final int highestValidCSPidx = maxIdxAvailCSPsols();
//        
//        int trials = 0;
//        boolean bMoved = false;
//        
//        for (int j = 0; j < this.tourSize_; j++) {
//            //directions = MyAlgorithms.getDirection(candidates.get(j).getValsCopy(), candidates.get((j+1)%tourSize_).getValsCopy());
//            //maxDist = MyMath.norm(candidates.get(j).getValsCopy(), candidates.get((j+1)%tourSize_).getValsCopy(), MyMath.DIST_EUCLEADIAN);
//        
//            
//            p1 = candidates.get(j);
//            
//            if(highestValidCSPidx >= 0 && Math.random()<0.5 && !p1.isSolution()){                        
//                pickedCons = MyRandom.randperm(0, p1.getRankComponents().size()-1).get(0);
//                pickedCons = p1.getRankComponents().get(pickedCons);
//                
//                p2 = new Chromosome(this.userInput_.solutionBy,  userInput_);
//                randPickIdx = MyRandom.randperm(0, highestValidCSPidx).get(0);
//                p2.setVals(CSPsols.get(randPickIdx).get(pickedCons), maxCSPval, CSPsols);   
//                
//            }else{
//                if(highestValidCSPidx == -1){ //I use then when it is difficult to find CSP
//                    for (int i = 0; i < p1.getRankComponents().size(); i++) {
//                        pickedCons = p1.getRankComponents().get(i);
//                        
//                        
//                        if(CSPsols.get(0).get(pickedCons).isEmpty()){
//                            p2 = null;
//                            continue;
//                        }else{
//                            p2 = new Chromosome(this.userInput_.solutionBy,  userInput_);
//                            p2.setVals(CSPsols.get(0).get(pickedCons), maxCSPval, CSPsols);
//                            break;
//                        }
//                    }
//                }                               
//                
//                if(p2 == null)
//                    p2 = (Chromosome)candidates.get((j+1)%tourSize_).clone();                
//            }            
//            
//
//            delta = MyMath.vectorSubtraction(p2.getValsCopy(), p1.getValsCopy());
//            int k = 0;
//            ArrayList<Double> prevVals = null, newVals=null;
//            int forceFind = 0;
//            do{
//                for (k = 1; k <= move; k++){ //) move; k++) {
//                    //find which direction to move?
//                    childChrome = new Chromosome(this.userInput_.solutionBy, this.userInput_);
//                    //approachDist = Math.pow(ratio,k)*maxDist;                
//                    newDelta = MyMath.constMultiplicationToVector(Math.pow(bringCloserRatio,k), delta); 
//                    newVals = MyMath.vectorAddition(p1.getValsCopy(), newDelta);
//                    childChrome.setVals(newVals, maxCSPval, CSPsols);
//
//                    //vp = new VirusProliferate(movingChrome.vals.toArray(), this.range_);
//
//                    //**************************************************************************************************************//
//                    //NOTE: ofsp changed <= sign to < sign
//                    //It is now giving me less solutions
//                    //It is good or bad......... I don't know.... It only promotes local search.
//                    //if(childChrome.getRank() < p1.getRank()){// || (childChrome.getRank() <= p1.getRank() && move == 1)){                    
//                        bMoved = false;
//                        
//                        if(!childChrome.myParent(p1)){ //there is a gap/black hole so no point searching further.
//                            bMoved = false;
//                            break;
//                        }
//                        if(childChrome.isSolution() || p1.isMyChild(childChrome)){
//                            offspring.add(childChrome); // if using do while then put this line in correct place.
//                            if(childChrome.isSolution()){
//                                ;//throw new SolutionFoundException("Sol found during crossover...");
//                            }
//                            bMoved = true;
//                            break;
//                        }
//                   // }
//
//                    prevVals = newVals;
//
//                }
//                if(bMoved && k>=1){
//                    p1  = childChrome;
//                    delta = MyMath.vectorSubtraction(p2.getValsCopy(), p1.getValsCopy());
////                }else if(bMoved && k>=2){
////                    p1.setVals(prevVals, maxCSPval, CSPsols);
////                    p2.setVals(newVals, maxCSPval, CSPsols);
////                    delta = MyMath.vectorSubtraction(p2.getValsCopy(), p1.getValsCopy());
//                }else{
//                    forceFind = 1;
//                }
//                forceFind++;
//            }while(forceFind < 1);
//            //chk boundary...
//
//            
////            if(!bMoved && trials < 5 && highestValidCSPidx >= 0){
////                trials++;
////                j--;
////            }
////            
////            if(!bMoved && trials >= 5 && highestValidCSPidx >= 0){
////                trials = 0;
////                candidates.set(++j, p2);
////            }
//            
//        }
//        return offspring;
//    }



    /**
     * Mutate the given set of offspring.
     * @param offspring mutation applied only to offspring
     */
    private void mutation(ArrayList<Chromosome> offspring) throws UnsupportedDataTypeException, SolutionFoundException{
        if(offspring.isEmpty()){
            return;
        }
        if(!userInput_.doMutation){
            return;
        }

         //update the offspring
        if(this.dataType_.contains("Integer")){
            mutationInteger(offspring);
        }
        else if(this.dataType_.contains("Double")){
            mutationDouble(offspring);
        }else{
            throw new UnsupportedDataTypeException("Only supports Integer and Double data type");
        }
    }

     /**
     * mutationDouble only mutate Doubles. It uses Polynomial Mutation as described in NSGA - II <br>
     * <B>Note</B> that offspring ArrayList is updated here.
     * @param offspring offspring generated after crossover.
     */
    private void mutationDouble(ArrayList<Chromosome> offspring) throws SolutionFoundException{
        int size = offspring.size();
        ArrayList<Integer> randInts;
        Chromosome temp;
        double val;
        double rand;
        double add;
        int muteBits = (int)Math.ceil(0.1*userInput_.totalDecisionVars); //10%
        
        for (int i = 0; i < size; i++) {
            
            try{
                if(Math.random()<1.0/userInput_.totalDecisionVars){
                    randInts = MyRandom.randperm(0, size-1);
                    temp = offspring.get(randInts.get(0));                   

                    for (int j : MyRandom.randperm(0, userInput_.totalDecisionVars-1).subList(0, muteBits)){                   
                    //for (int j = 0; j < userInput_.totalDecisionVars; j++) {
                        val = temp.getVals(j);
                        rand = Math.random();
                        if(rand<0.5)
                            add = Math.pow(2.0*rand,1.0/(MUM+1)) -1;
                        else
                            add = 1- Math.pow(2.0*(1-rand),1.0/(MUM+1));

                        val = val+add;

                        if(val>userInput_.maxVals.get(j))
                            val = userInput_.maxVals.get(j);
                        else if(val<userInput_.minVals.get(j))
                            val = userInput_.minVals.get(j);

                        temp.replaceVal(j, val, maxCSPval); 
                    }               
                }  
            }catch (SolutionFoundException sfe){
                throw sfe;
            }catch(Exception e){
                e.printStackTrace();
            }
        }
    }
    
    
    private void mutationSwap(Chromosome ch, double maxSwap, final double maxIteration) throws SolutionFoundException{
        ArrayList<Integer> randVal; // = MyRandom.randperm(0, ch.getValsCopy().size()-1);
        ArrayList<Double> vals = new ArrayList<Double>();
        double val0, val1;
        int lowIdx, hiIdx;
        //double maxSwap = 1; //0.05*ch.getValsCopy().size();
//        final double maxIteration = 4;
        double bfFitness, afFitness;
        
        if(ch.getValsCopy().size()<2){
            return;
        }
        
        bfFitness = ch.getFitnessVal(0);
        
        for (int j = 0; j < maxIteration; j++) {
            for (int i = 0; i < maxSwap; i++) {
                randVal = MyRandom.randperm(0, ch.getValsCopy().size()-1);

                lowIdx = randVal.get(0);
                hiIdx = randVal.get(1);

                if(lowIdx > hiIdx){
                    lowIdx = randVal.get(1);
                    hiIdx = randVal.get(0);
                }

                val0 = ch.getVals(lowIdx);
                val1 = ch.getVals(hiIdx);

                vals.add(val0);
                vals.add(val1);

                try {
                    ch.remove(lowIdx);
                    ch.remove(hiIdx-1);
                    //ch.appendVal(val1);
                    //ch.appendVal(val0);
                } catch (Exception e) {
                    e.printStackTrace();    
                    System.out.println("uffff... yeh ulfat..");

                }   
            }

            for (Double v : vals) {
                ch.appendVal(v, maxCSPval);
            }
        
            afFitness = ch.getFitnessVal(0);
            
            if(afFitness!=bfFitness){
                if(afFitness<bfFitness){
                    afFitness = afFitness;
                }
                break;
            }
        }
    }
    
    private void mutationSwapNew(Chromosome ch, double maxSwap){
        ArrayList<Integer> randVal; // = MyRandom.randperm(0, ch.getValsCopy().size()-1);
        ArrayList<Double> vals = new ArrayList<Double>();
        double val0;
        int lowIdx;
        //double maxSwap = 1; //0.05*ch.getValsCopy().size();
//        final double maxIteration = 4;
        
        if(ch.getValsCopy().size()<2){
            return;
        }
        
        
        for (int i = 0; i < maxSwap; i++) {
            randVal = MyRandom.randperm(0, ch.getValsCopy().size()-1);

            lowIdx = randVal.get(0);               

            val0 = ch.getVals(lowIdx);
            vals.add(val0);

            try {
                ch.remove(lowIdx);
            } catch (Exception e) {
                e.printStackTrace();    
                System.out.println("uffff... yeh ulfat..");
            }   
        }

            
        ch.tryForcedCSPsolUpdate();
       
    }
    
    private void mutationGroupSwap(Chromosome ch) throws SolutionFoundException{
        ArrayList<Integer> randVal = MyRandom.randperm(0, ch.getSatisfaction().size()-1);
        int Idx0, Idx1;
        int loc;
        ArrayList<ArrayList<Integer>> temp = new ArrayList<ArrayList<Integer>>();
        final int sz = ch.getSatisfaction().size();
        
        if(Math.random()<0.1){
            loc = MyRandom.randperm(0, sz -1).get(0);
            
            for (int i = 0; i < sz; i++) {
                temp.add(ch.getSatisfaction().get(i));
            }
            
            for (int i = 0; i < sz; i++) {
                ch.getSatisfaction().set(i,temp.get((i+loc)%sz));
            }
           
        }else{
            Idx0 = randVal.get(0);
            Idx1 = randVal.get(1);

            ArrayList<Double> list0 = ch.getSatisfaction().get(Idx0);
            ArrayList<Double> list1 = ch.getSatisfaction().get(Idx1);

            ch.getSatisfaction().set(Idx0, list1);
            ch.getSatisfaction().set(Idx1, list0);
        }

        //ch.refreshValVsConstIdx();       
        ch.refreshFitness(maxCSPval);        
    }
    
    /**
     * METHOD IS NOT TESTED. TEST IT FIRST BEFORE USE.
     * mutationInteger only mutate integers. It uses swap elements technique. 
     * so that it disrupts order more to get new allel values<br>
     * <B>Note:</B> that offspring ArrayList is updated here but rank will remain
     * same because swapping satisfaction value will produce same result. It may
     * only give different results in crossover.
     * @param offspring offspring generated after crossover.
     */
    private void mutationInteger(ArrayList<Chromosome> offspring) throws SolutionFoundException{
        ArrayList<Integer> randDim;
        ArrayList<Integer> randVal;
        Double temp = 0.0;
        int muteBits;
        
        //System.out.println("testing... " + offspring);

        if(userInput_.domainVals == null || userInput_.domainVals.isEmpty()){ //mutation not supported
            return;
        }

        if(externalData_ == null){ //currently works only for external data
            return;
        }
        
        //Technique 1: swapping values
        //<<
//        for (int ofsp = 0; ofsp < offspring.size(); ofsp++) {            
//            if(Math.random()<1.0/offspring.get(ofsp).getValsCopy().size()){
//                //Only deal with valid values...
//                randDim = MyRandom.randperm(0, offspring.get(ofsp).getRankComponents().size()-1);                
//                
//                if(randDim.size()<2){ //swapping not possible
//                    continue;
//                }else{
//                    holdVals = (ArrayList<Double>)offspring.get(ofsp).getValsCopy().clone();
//                    holdVals.set(randDim.get(0), offspring.get(ofsp).getValsCopy(randDim.get(1)));
//                    holdVals.set(randDim.get(1), offspring.get(ofsp).getValsCopy(randDim.get(0)));                    
//                    offspring.get(ofsp).setVals(holdVals);                    
//                }
//            }
//        }
        //>>
        
        
        ArrayList<Double> vals;
        //ArrayList<Double> noGoods; 
        int expectedVal;
        ArrayList<Double> noGoods;
        
        
        
        //Technique 2: mutate a given value from available domain value;
        //<<
//        for (int ofsp = 0; ofsp < offspring.size(); ofsp++) {  
        for (Chromosome offsp : offspring) {                    
            if(Math.random()<MUTATION_RATE){ //1.0/userInput_.totalDecisionVars){ //>1.0/offspring.get(ofsp).getValsCopy().size() || bStagnant){
                               
                vals = offsp.getValsCopy();
                Collections.sort(vals);
                expectedVal = 0;
//                noGoods = new ArrayList<Double>();
//
//                for (int i = 0; i < vals.size(); i++) {                        
//                    if(vals.get(i).intValue() != expectedVal){
//                        for (int j =expectedVal; j < vals.get(i).intValue(); j++) {
//                            noGoods.add(j*1.0); 
//                        } 
//                        expectedVal = vals.get(i).intValue();
//                    } 
//                    expectedVal++;
//                }
//
//                for (int i = vals.get(vals.size()-1).intValue()+1; i < userInput_.totalDecisionVars; i++) {
//                    noGoods.add(i*1.0); 
//                }
                vals.clear();
                noGoods = offsp.findNoGoods();
                if(noGoods.isEmpty()){
                    continue; //nothing to replace with
                }
                
                if(offsp.getValsCopy().size()<2){ //swapping not possible
                    continue;
                }else{
                    if(userInput_.domainVals == null){
                        continue;
                    }else if(userInput_.domainVals.isEmpty()){
                        continue;
                    }
                    
                    muteBits = 1;
                    if(bStagnant){
                        muteBits = Math.max(1,(int)(offsp.getValsCopy().size()*0.2));
                    }
                    //randVal = MyRandom.randperm(0,offsp.noGoods.size()-1);
                    
                    for (int j = 0; j < muteBits && j<noGoods.size(); j++) {
                        
                        randVal = MyRandom.randperm(0,noGoods.size()-1);
                        
                        if(bStagnant){ //Important... must refresh in every iteration....
                            muteBits = Math.max(1,(int)(offsp.getValsCopy().size()*0.2));
                        }

                        //Only deal with valid values...
                        randDim = MyRandom.randperm(0, offsp.getValsCopy().size()-1);
                        
                        if(randDim.get(0) >= offsp.getValsCopy().size()){
                            System.out.println("ee kaisey sake...");
                            System.out.println(vals);
                        }
                        try{
//                            randVal = MyRandom.randperm(0,userInput_.domainVals.get(randDim.get(0)).size()-1);
//                            randVal = MyRandom.randperm(0,noGoods.size()-1);
                            if(!externalData_.isHighlyConstrained(offsp.getVals(randDim.get(0)).intValue())) //in optimization mode noGoods is empty so automatically this won't be executed.
                                offsp.replaceVal(randDim.get(0),noGoods.get(randVal.get(0)),maxCSPval);

//                            int prevVal =  vals.get(randDim.get(0)).intValue();
//                            for (int k = 0; k < userInput_.domainVals.get(prevVal).size(); k++) {
//                                temp = userInput_.domainVals.get(prevVal).get(randVal.get(k));
//                                if(vals.get(randDim.get(0)) != temp){                                
//                                    offspring.get(ofsp).replaceVal(randDim.get(0), temp);//Be warned! may create duplicate values..
//                                    break;
//                                }   
//                            }
                        
                        }catch(Exception e){
                            e.printStackTrace();
                            System.out.println("arey??");
                        }
                        
                    }
                             
                }
            }
        }
        //>>
    }

    private void sortAndReplace(int gen) throws Exception, SolutionFoundException{
        if (userInput_.dataType.contains("Integer")){
            noViolationSortAndReplaceInteger(gen); //duplicateSatisfactionSortAndReplace();
        }else if (userInput_.dataType.contains("Double")){
            throw new Exception("Not supported in this version...");
        }
        else{    
            throw new Exception("Incorrect use of data types");
        }
    }
 
    
    private void noViolationSortAndReplaceInteger(int gen) throws SolutionFoundException, Exception{        
        ArrayList<Chromosome> sols = new ArrayList<Chromosome>();
        ArrayList<Chromosome> nonSols = new ArrayList<Chromosome>();               
        final int funcionalConstraints = userInput_.total__updatedConstraints; //userInput_.totalConstraints - userInput_.totalDecisionVars+1;
        ArrayList<ArrayList<Chromosome>> grouping = new ArrayList<ArrayList<Chromosome>>();
        int front2[] = new int[2];
        
        int slotSize;
        int FirstSlotSize;
        int empty;
        int incompleteSlots;
        int slotAddition;
        int satis;
        
//        for (Chromosome ch : chromosomes_) {
//            ch.refreshFitness(maxCSPval);
//        }
        
//        curBest_ = bestSoFar.getRankComponents().size();
        
        
        
        if(bOptimizationMode)
            curBest_ = MyMath.roundN(bestSoFar.getFitnessVal(0),FIT_DP_LIMIT);
        else
            curBest_ = MyMath.roundN(bestSoFar.getFitnessVal(0),FIT_DP_LIMIT);
        
        if(prevBest_ == curBest_){
            stillSameBestCount++;
        }else{
            stillSameBestCount = 0;
        }          
        
        if(stagnantVisit >= 3 || stillSameBestCount == 0){
            stillSameBestCount = 0;
            bStagnant = false;
            stagnantVisit = 0;
        }        
        
        Chromosome ch; 
        for (int i = 0; i < chromosomes_.size(); i++) {
            ch = chromosomes_.get(i);
            
//            if(Chromosome.tmpSortBy != userInput_.solutionBy){
//                System.out.println("bad robot...");
//            }
            
            if(ch.isSolution()){
                sols.add(ch);
                chromosomes_.remove(i);
                i--;
            }else{
                //ch.tempSortBy = Chromosome.BY_VIOLATIONS; //sortBy; //
                nonSols.add(ch);
                chromosomes_.remove(i);
                i--;
            }
        }

        chromosomes_.clear();  
        
        int solPop = Math.min(userInput_.population/2,sols.size());
        int nonSolPop = Math.min(userInput_.population/2,nonSols.size()); //userInput_.population - solPop;
        
        if(solPop < userInput_.population/2){
            nonSolPop = userInput_.population - solPop;
        }
        
        if(nonSolPop < userInput_.population/2){
            solPop = userInput_.population - nonSolPop;
        }
                
        
        if(!sols.isEmpty()){
            bOptimizationMode = true;         
            this.maxCSPval = sols.get(sols.size()-1).getFitnessVal(1);                                           
        }
        
        if(bOptimizationMode){
            this.userInput_.bWeighted = true; 
            for (Chromosome newc : sols) {                        
                newc.refreshFitness(maxCSPval);
            }
            Chromosome.tmpSortBy = Chromosome.BY_FITNESS; //userInput_.solutionBy;
            Collections.sort(sols);
            bestSoFar = (Chromosome)sols.get(0).clone();
            
            this.userInput_.bWeighted = false; 
            for (Chromosome newc : nonSols) {                        
                newc.refreshFitness(maxCSPval);
            }
            this.userInput_.bWeighted = true; 
            
            for (Chromosome rep : suspended_) {                        
                rep.refreshFitness(maxCSPval);
            }
        }         

        try{ 
            if(!nonSols.isEmpty()){   

                
                for (Chromosome nsch : nonSols) {
                    if(nsch.getFitnessVal(0) - userInput_.total__updatedConstraints == 0){ //checking in nonsols means it is certain these are transitions not the final sols
                        break;
                    }
                }
                
                //<<
                nonSols = fitnessThenNoveltySort(nonSols, nonSols.size());
                double ratio = 0.1;
                if(userInput_.totalConstraints<10){
                    ratio = 0.5;
                }
                final int slots = Math.min((int)(ratio*userInput_.totalConstraints), nonSolPop); //0.5 or 1.0 for small ones otherwise 0.1
                final int remainder = (int)Math.ceil(userInput_.totalConstraints/slots);
                slotSize = nonSolPop/(slots); //last one for infeasibles as well
                FirstSlotSize = slotSize + (nonSolPop - slotSize*(slots));
                empty = 0;
                incompleteSlots = 0;

                for (int i = 0; i < slots; i++) {//funcionalConstraints
                    grouping.add(new ArrayList<Chromosome>());
                }

                satis = -1;                
//                    for (Integer i : MyRandom.randperm(0, nonSols.size()-1)) {
                for(int i = 0; i<nonSols.size(); i++){
                    satis = Math.min(funcionalConstraints-1,nonSols.get(i).getFitnessVal(0).intValue()); //1 to total constraints
                    satis = Math.min((int)Math.floor(satis/remainder),slots-1);
                    grouping.get(satis).add(nonSols.get(i));
                }

                Chromosome.tmpSortBy = userInput_.solutionBy;//Chromosome.BY_RO;//userInput_.solutionBy;//
                for (ArrayList<Chromosome> grp : grouping) {
                    Collections.sort(grp);
                }

                int fr = 0;                                        
                for (int i = 0; i < grouping.size(); i++) {
                    ArrayList<Chromosome> g = grouping.get(i); 

                    if(g.size()>0 && fr < 2){
                        front2[fr++] = i;
                    }                        
                }

                empty = 0;
                incompleteSlots = 0;
                for (int i = 0; i < grouping.size(); i++) {
                    if(grouping.get(i).size()<slotSize){
                        empty += slotSize-grouping.get(i).size();
                        incompleteSlots++;
                    }
                }

                slotAddition = empty/(slots-incompleteSlots); //empty slot space has to be distributed to filled/partially filled slots 
                slotSize += slotAddition;
                FirstSlotSize += slotAddition;
                FirstSlotSize += empty - (slots-incompleteSlots)*slotAddition;

                ArrayList<Chromosome> additionals = new ArrayList<Chromosome>();

                boolean bFirstSlotAdded = false;
                for (int i = 0; i < grouping.size(); i++) {
                    if(!bFirstSlotAdded && grouping.get(i).size() >= FirstSlotSize){

                        if(grouping.get(i).size() > FirstSlotSize )
                            additionals.addAll(grouping.get(i).subList(FirstSlotSize, grouping.get(i).size()));

                        grouping.set(i, new ArrayList<Chromosome>(grouping.get(i).subList(0, FirstSlotSize)));

                        bFirstSlotAdded = true;
                        continue;
                    }

                    if(grouping.get(i).size() > slotSize )
                        additionals.addAll(grouping.get(i).subList(slotSize, grouping.get(i).size()));

                    grouping.set(i, new ArrayList<Chromosome>(
                            grouping.get(i).subList(
                                0, Math.min(
                                    slotSize,grouping.get(i).size()))));
                }

                int sz = 0;
                for (ArrayList<Chromosome> g : grouping) {
                    sz += g.size();                    
                }

                if(additionals.size() + sz != nonSols.size()){
                    sz =sz;
                }

                nonSols.clear();

                //int fr = 0;                                        
                for (int i = 0; i < grouping.size(); i++) {
                    ArrayList<Chromosome> g = grouping.get(i); 
                    nonSols.addAll(g);                        
                }

                nonSols.addAll(additionals.subList(0, nonSolPop-nonSols.size())); 
                Chromosome.tmpSortBy = Chromosome.BY_FITNESS;
                
                
                //>>

                
                 if((stillSameBestCount >= SAME_BEST_GENERATIONS)){
                    bStagnant = true;                   
                    stagnantVisit++;           
                    
                    if(stagnantVisit == 1){
                        System.out.println("*** temp -- removed...");
                        System.out.println("cur best: " + bestSoFar.getSatisfaction());
//                        ArrayList<Double> nearestNeighbor = null;
//                        tabuDist = -1.0;
//                        satis = nonSols.get(0).getRankComponents().size();                                                
//                        if(grouping.size()>=2){
//                            if(grouping.get(front2[0]).size()>0 && grouping.get(front2[1]).size()>0){ //obviously... huun..
//                                tabuDist = Math.abs(MyMath.norm(grouping.get(front2[0]).get(0).getValsCopy(), 
//                                        grouping.get(front2[1]).get(0).getValsCopy(),MyMath.DIST_EUCLEADIAN));
//                            }
//                        }                        
//                                             
//                        if(tabuDist>0){
                            addDynamicConstraint(nonSols.get(0).getSatisfaction(), tabuDist);                        
//                        }
                        for (Chromosome c : nonSols) {
                            c.refreshFitness(maxCSPval);
                        }  
                        
                        Chromosome.tmpSortBy = Chromosome.BY_FITNESS;
                        Collections.sort(nonSols);
                    }                            
                }
                
                if(!bOptimizationMode){
                    bestSoFar = (Chromosome)nonSols.get(0).clone(); //It is not necessary that get(0) it the best accroding to fitness val                                                                    // but it donsen't matter in case of CSP...
                }

            }
        } catch(Exception e){
            throw e;
        }    
        
        sols = new ArrayList<Chromosome>(sols.subList(0, solPop));
        nonSols = new ArrayList<Chromosome>(nonSols.subList(0, nonSolPop));

        chromosomes_.clear(); //just to be safe;
        chromosomes_.addAll(sols);
        chromosomes_.addAll(nonSols);
            
        if(chromosomes_.size() != userInput_.population){
            System.err.println("population size error on noViolationSortAndReplace.");
            Application.getInstance().exit();
        }

        
        Chromosome.tmpSortBy = Chromosome.BY_FITNESS;//userInput_.solutionBy;
        
 
        prevBest_ = curBest_;
        if(bStagnant){
            randomDeath(gen, 1, true); //spacre top one            
        }
        else            
            randomDeath(gen, (int)(userInput_.population*0.1),false);//,true  
        

        //<<incrementalty...
        if(bestSoFar.getFitnessVal(0) == 0){
            if(bestSoFar.isSolution()){
                throw new SolutionFoundException("solution found: not properly updated");
            }
            
            dynamicConstraintNo++;            
            if(dynamicConstraintNo>= MAX_FUNCTIONAL_CONSTRAINTS){
                negFeasibleRange++;
                dynamicConstraintNo = MAX_FUNCTIONAL_CONSTRAINTS;
            }
                        
            for (Chromosome c : chromosomes_) {
                c.refreshFitness(maxCSPval);
            }
            bOptimizationMode = false;
            this.userInput_.bWeighted = false;
            stillSameBestCount = 0;
            bStagnant = false;
            stagnantVisit = 0;

            Chromosome.tmpSortBy = Chromosome.BY_FITNESS;//userInput_.solutionBy;
            Collections.sort(chromosomes_);
            bestSoFar = (Chromosome)chromosomes_.get(0).clone();            
        }
        //>>..........
        //??? not sure if this is required
//        for (Chromosome chr : chromosomes_) {
//            chr.refreshFitness(maxCSPval);
//        }
        
//        if(gen>5)
//        for (Chromosome chr : chromosomes_) {
//            chr.refreshFitness(maxCSPval);
//            chr.setVals(new ArrayList<Double>(chr.getValsCopy().subList(Math.min(1, chr.getValsCopy().size()), chr.getValsCopy().size())), maxCSPval);
//            chr = chr;
//            //
//        }
    }
    
    /**
     * NOTE: 
     * @param center
     * @param radius 
     */
    private void addDynamicConstraint(ArrayList<ArrayList> center, double radius){ 
        //negFeasibleRange = 0; no need 
        //dynamicConstraintNo = 0;
        //userInput_.totalConstraints++; done inside externalData_.addTabuConstraint
        //MAX_FUNCTIONAL_CONSTRAINTS = userInput_.totalConstraints - userInput_.totalDecisionVars; don't need here
        
        if(externalData_ != null)
            externalData_.addTabuConstraint(center); //center = local optimal position 
//        else
//            CCSPfns.addTabuConstraint(center, radius);    
        //MUST...
        for (Chromosome c : chromosomes_) {
            c.cleanFitnessHistory();
        }
    }
    
    /**
     * NOTE: This function <B>DOES NOT</B> sort the chromosomes inside a ranked group.
     * If it has <I>n</I> ranks/groups, it only tries to give best ranked chromosomes,
     * then the leftovers are <B>ONLY</B> sorted according to fitness.
     * You must sort each ranked groups separately afterwards.
     * @param in
     * @param size
     * @return 
     */
    private ArrayList<Chromosome> fitnessThenNoveltySort(final ArrayList<Chromosome> in, final int size){        
        Chromosome.tmpSortBy = Chromosome.BY_FITNESS;
        Collections.sort(in);
        
        ArrayList<Chromosome> out = null;
        final int minAcceptedRank = in.get(size-1).getFitnessVal(0).intValue(); //getRankComponents().size(); //violations
         
        int safePointer = -1;
        Chromosome chrome;
        ArrayList<Chromosome> temp = new ArrayList<Chromosome>();
        
        try{
            if(in.size()<=1 || in.size() <= size){
                out = in;
                throw new ExecutionException(null);
            }

            for (Chromosome chrm : in) {
                if(chrm.getFitnessVal(0) == minAcceptedRank){
                    chrome = (Chromosome)chrm;
                    temp.add(chrome); //exract chromosomes with max accepted violations.
                }else if(chrm.getFitnessVal(0) > minAcceptedRank){
                    safePointer++;
                }else{
                    break;
                }
            }

            if(size <= safePointer+1){
                out = new ArrayList<Chromosome>(in.subList(0, size));//not get only required sorted ones.  
            }else{
                out = new ArrayList<Chromosome>(in.subList(0, safePointer+1));//not get only required sorted ones.                        
                Chromosome.tmpSortBy = Chromosome.BY_RO; //userInput_.solutionBy;
                if(!temp.isEmpty()){
                    temp = temp;
                }
                Collections.sort(temp);
                out.addAll(temp.subList(0, size-safePointer-1));                
                out = new ArrayList<Chromosome>(out.subList(0, size));
            }        
        }catch(ExecutionException ee){
            out = out;
        }
        catch(Exception e){
            e.printStackTrace();
            out = out;
        }

        Chromosome.tmpSortBy = userInput_.solutionBy;
        return out;
    }
    
    
    private void randomDeath(int gen, int spareSize, boolean bImportSuspended) throws SolutionFoundException{
        int d;
        ArrayList<Chromosome> newRandPop = new ArrayList<Chromosome>();
        d = (int)Math.round(this.REPLACE_PERCENT*userInput_.population);
        
//        try {
//            initializeChromosomes(newRandPop, d, gen);             
//        } catch (Exception e) {
//            e.printStackTrace();
//            Application.getInstance().exit();
//        }   


//            hey??????????? its not a random............. bluffffffffff
//            for (int ofsp = 0; ofsp < d; ofsp++) {  //spare the top d chromes          
//                chromosomes_.set(userInput_.population-1-ofsp, newRandPop.get(ofsp));
//            }
            
            
//            int temp = 0;
//            Chromosome ch;
//            for(int i: MyRandom.randperm(spareSize, chromosomes_.size()-1).subList(0, d)){
//                if(suspended_.size()<userInput_.population && chromosomes_.get(i).isFeasible() && !bStagnant)
//                    suspended_.add(chromosomes_.get(i)); //getting reference? it is ok as it will be deleted below.
//                
//                if(gen<5 || !bImportSuspended || suspended_.isEmpty()){//0.01*userInput_.generation){
//                    chromosomes_.set(i,newRandPop.get(temp));
//                    temp++;
//                }
//                else{    
//                    ch = suspended_.remove();
//                    //ch.refreshFitness(maxCSPval);
//                    chromosomes_.set(i,ch);   
//                    suspended_.add((Chromosome)ch.clone());
//                }
//            }        

    }
    
    
//    private void noViolationSortAndReplace(int gen) throws Exception{
//        ArrayList<Chromosome> diverse = new ArrayList<Chromosome>();
//        ArrayList<Chromosome> newRandPop = new ArrayList<Chromosome>();
//        Chromosome chrome;
//        double maxAcceptedRank;
//        int safePointer = -1;
//        int d;
//
//
//        for (Chromosome chrm : chromosomes_) {
//            chrm.tempSortBy = Chromosome.BY_FITNESS;
//        }
//        
//        Collections.sort(chromosomes_);//sorted according to violation
//        
//        bestSoFar = chromosomes_.get(0);
//        
//        for (Chromosome chrm : chromosomes_) {
//            if(chrm.isSolution()){
//                throw new SolutionFoundException("All constraints satisfied");
//            }
//        }
//
//
//        int tempSize;
//        ArrayList<Double> tempVals;
//
//        bStagnant = false;
//        if(gen>userInput_.generation*0.05 && gen%10 == 0){
//            bStagnant = true;
//            System.diverse.println("*** diverse -- removed...");
//            tempSize = (int)(1*SAME_BEST_VAL_PERCENT*userInput_.population);        
//            
//            for (int ofsp= 0; ofsp < tempSize; ofsp++) {
//                tempVals = chromosomes_.get(ofsp).negateVals();
//                chromosomes_.get(ofsp).setVals(tempVals);
//            }            
//            
//            Collections.sort(chromosomes_);          
//        }
//
//        maxAcceptedRank = chromosomes_.get(userInput_.population-1).getFitnessVal(0);
//
//        try{                       
//            for (Chromosome chrm : chromosomes_) {
//                if(MyMath.roundN(chrm.getFitnessVal(0),2) == MyMath.roundN(maxAcceptedRank,2)){
//                    chrome = (Chromosome)chrm;
//                    chrome.tempSortBy = userInput_.solutionBy ;
//                    chrome.tempRo = getRoValue(chrm);
//                    diverse.add(chrome); //exract chromosomes with max accepted violations.
//                }else if(MyMath.roundN(chrm.getFitnessVal(0),2) < MyMath.roundN(maxAcceptedRank,2)){
//                    safePointer++;
//                }else{
//                    break;
//                }
//            }
//            chromosomes_ = new ArrayList<Chromosome>(chromosomes_.subList(0, safePointer+1));//not get only required sorted ones.
//
//            for (Chromosome chrm : diverse) {
//                chrm.tempSortBy = Chromosome.BY_RO; //MUST DO before sorting
//            }
//            Collections.sort(diverse);
//            for (Chromosome chrm : diverse) {
//                chrm.tempSortBy = Chromosome.BY_FITNESS; //MUST DO before sorting
//            }
//
//            chromosomes_.addAll(diverse.subList(0, userInput_.population-safePointer-1));
//
//            if(chromosomes_.size() != userInput_.population){
//                System.err.println("population size error on noViolationSortAndReplace.");
//                Application.getInstance().exit();
//            }
//            
//            for (Chromosome c : chromosomes_) {
//                c.tempSortBy = userInput_.solutionBy;
//            }
//
//            randomDeath();
//        }catch(MyException me){
//            me.showMessageBox();
//        }
//        catch(Exception e){
//            throw e;
//        }
//    }        
    
    private boolean XXXforceFindSolution(Chromosome chrom){                
        ArrayList<Double> vals = chrom.getValsCopy();
        ArrayList<Double> noGoods = new ArrayList<Double>();
       
        Collections.sort(vals);
        int expectedVal = 0;
        
        
        //WRONG CODE..... SEE tryForcedCSPsolUpdate IN TT
        //<<
        for (int i = 0; i < vals.size(); i++) {                        
            if(vals.get(i).intValue() != expectedVal){
                for (int j =expectedVal; j < vals.get(i).intValue(); j++) {
                    noGoods.add(j*1.0); 
                } 
                expectedVal = vals.get(i).intValue();
            } 
            expectedVal++;
        }
        
        for (int i = vals.get(vals.size()-1).intValue()+1; i < userInput_.totalDecisionVars; i++) {
            noGoods.add(i*1.0); 
        }
        //>>

        int prevLength;
        int newLenght=-1;
        Double removedVal;
        
        try{                    
            for (Double ng : noGoods) {               
                prevLength = chrom.getValsCopy().size();
                for (int i = 0; i < chrom.getValsCopy().size(); i++) {
                    removedVal = chrom.getVals(i);
                    if(externalData_.isViolated(removedVal.intValue(), ng.intValue()) && !externalData_.isHighlyConstrained(removedVal.intValue())){
                        chrom.replaceVal(i, ng, maxCSPval);
                        chrom.appendVal(removedVal, maxCSPval);
                    }
                    newLenght = chrom.getValsCopy().size();

                    if(newLenght>prevLength){//successful addition
                        break;
                    }
                    if(newLenght<prevLength){//successful addition
                        System.out.println("keee kaisey sake?????");
                    }
                    
                }            
            }
        }catch(Exception e){
            e.printStackTrace();
        }
        
        if(newLenght == userInput_.totalDecisionVars){
            return true;
        }else{
            return false;
        }
        
    }
    
    
//    private void xxx(int gen) throws Exception{
//        for (Chromosome chrm : chromosomes_) {
//            chrm.tempRo = getRoValue(chrm);
//            chrm.tempSortBy = userInput_.solutionBy ;//Chromosome.BY_VIOLATIONS; //MUST DO before sorting
//        }
//        
//        Collections.sort(chromosomes_);//sorted according to violation/satisfactions
//        
//        bestSoFar = chromosomes_.get(0);
//        curBest_ = bestSoFar.getRank();
//        
//        for (Chromosome chrm : chromosomes_) {        
//            if(chrm.isSolution()){
//                bestSoFar = chrm;
//                bestSoFar.tempSortBy = userInput_.solutionBy;
//                throw new SolutionFoundException("All constraints satisfied");
//            }
//            chrm.tempSortBy = Chromosome.BY_RO;    
//        }
//
//        Collections.sort(chromosomes_);  
//        chromosomes_ = new ArrayList<Chromosome>(chromosomes_.subList(0, userInput_.population));
//        
//        for (Chromosome chrm : chromosomes_) {
//            chrm.tempSortBy = userInput_.solutionBy ;//Chromosome.BY_VIOLATIONS; //MUST DO before sorting
//        }
//    }
    
    
//    private void setUniqueChromosomes(){
////        for (int i = 0; i < chromosomes_.size()-1; i++) {
////            for (int j = i+1; j < chromosomes_.size(); j++) {
////                if(chromosomes_.get(i).equals(chromosomes_.get(j))){
////                    chromosomes_.remove(j);
////                    j--;
////                }
////            }
////        }
//        
//        Collections.sort(chromosomes_);
//        
////        for (Chromosome chrm : chromosomes_) {
////            System.out.print(chrm.fitness_.get(0)+", ");
////        }
//        
//        for (int i = 0; i < chromosomes_.size()-1; i++) {
//            if(MyMath.roundN(chromosomes_.get(i+1).fitness_.get(0),2)
//                    == MyMath.roundN(chromosomes_.get(i).fitness_.get(0),2)){
//                chromosomes_.remove(i+1);
//                i--;
//            }            
//        }
//        
////        System.out.println("\n\n");
////        for (Chromosome chrm : chromosomes_) {
////            System.out.print(chrm.fitness_.get(0)+", ");
////        }
////        System.out.println("\n\n");
//    }
    
    
  
    
// <editor-fold defaultstate="collapsed" desc="Old commented code. May be useful :)">    
    
//////    /**
//////     * Sort based on violation preference. select the best ones and then use
//////     * ro values if same violation is found.
//////     * @throws Exception
//////     */
//////    private void noViolationSortAndReplace(int gen) throws Exception, SolutionFoundException{
//////        ArrayList<Chromosome> diverse = new ArrayList<Chromosome>();
//////        ArrayList<Chromosome> newRandPop = new ArrayList<Chromosome>();
//////        Chromosome chrome;
//////        int maxAcceptedViolation;
//////        int safePointer = -1;
//////        int d;
//////
//////        //Get immunity chromosomes...
//////        
//////        
////////        ArrayList<Chromosome> goodImmuneChrom = new ArrayList<Chromosome>();        
////////        for (Chromosome chrm : chromosomes_) {
////////            chrm.sortBy = Chromosome.BY_IMMUNITY;        
////////        }
////////        Collections.sort(chromosomes_);        
////////        goodImmuneChrom = new ArrayList<Chromosome>(chromosomes_.subList(0, (int)Math.round(this.IMMUNITY_PERCENT*userInput_.population)));
////////        chromosomes_.removeAll(goodImmuneChrom);
//////        
//////        for (Chromosome chrm : chromosomes_) {
//////            chrm.sortBy = userInput_.solutionBy ;//Chromosome.BY_VIOLATIONS; //MUST DO before sorting 
//////            if(chrm.isSolution()){
//////                throw new SolutionFoundException("All constraints satisfied");
//////            }
//////        }
//////                        
//////        Collections.sort(chromosomes_);//sorted according to violation
//////        bestSoFar = chromosomes_.get(0);
//////
//////        int immuneCount = 0;
////////        Chromosome tempChrome;
////////        
////////        chromosomes_.addAll(0, goodImmuneChrom);
////////        immuneCount = goodImmuneChrom.size();
////////        for (Chromosome chrm : chromosomes_) {
////////            chrm.sortBy = userInput_.solutionBy ;//Chromosome.BY_VIOLATIONS; //MUST DO before sorting 
////////        }
////////           
////////        for (Chromosome chrm : chromosomes_) {
////////            if(chrm.isSolution()){
////////                throw new SolutionFoundException("All constraints satisfied");
////////            }
////////        }
////////        
////////        for (int ofsp = 0; ofsp < chromosomes_.size(); ofsp++) {
////////            if(chromosomes_.get(ofsp).getImmunity() > 0){
////////                chromosomes_.get(ofsp).useImmunity();
////////                tempChrome = chromosomes_.remove(ofsp);
////////                chromosomes_.add(0, tempChrome);
////////                immuneCount++;
////////            }
////////            
////////            if(immuneCount>= (int)Math.round(this.IMMUNITY_PERCENT*userInput_.population))
////////                break;            
////////        }
//////    
//////        int tempSize;
//////        ArrayList<Double> tempVals;
//////
//////        bStagnant = false;
//////        if(gen>userInput_.generation*0.05 && gen%10 == 0){
//////            bStagnant = true;
//////            System.diverse.println("*** diverse -- removed...");
//////            tempSize = (int)(1*SAME_BEST_VAL_PERCENT*userInput_.population);
//////
//////            for (int ofsp= 0; ofsp < tempSize; ofsp++) {
//////                tempVals = chromosomes_.get(ofsp).negateVals();
//////                chromosomes_.get(ofsp).setVals(tempVals);
//////            }
//////            Collections.sort(chromosomes_);        
//////        }
//////
//////        maxAcceptedViolation = chromosomes_.get(userInput_.population-1).getRank();
//////
//////        try{  
//////            safePointer = immuneCount-1;
//////            for (int ofsp = immuneCount; ofsp < chromosomes_.size(); ofsp++) {
//////                if(chromosomes_.get(ofsp).getRank()==maxAcceptedViolation){
//////                    chrome = (Chromosome)chromosomes_.get(ofsp);
//////                    chrome.tempRo = getRoValue(chromosomes_.get(ofsp));
//////                    diverse.add(chrome); //exract chromosomes with max accepted violations.
//////                }else if(chromosomes_.get(ofsp).getRank() < maxAcceptedViolation){
//////                    safePointer++;
//////                }else{
//////                    break;
//////                }
//////            }            
//////            
//////            chromosomes_ = new ArrayList<Chromosome>(chromosomes_.subList(0, safePointer+1));//not get only required sorted ones.
//////
//////            for (Chromosome chrm : diverse) {
//////                chrm.sortBy = Chromosome.BY_RO; //MUST DO before sorting
//////            }
//////            Collections.sort(diverse);
//////            for (Chromosome chrm : diverse) {
//////                chrm.sortBy = userInput_.solutionBy; //MUST DO before sorting
//////            }
//////
//////            chromosomes_.addAll(diverse.subList(0, userInput_.population-safePointer-1));
//////
//////            if(chromosomes_.size() != userInput_.population){
//////                System.err.println("population size error on noViolationSortAndReplace.");
//////                Application.getInstance().exit();
//////            }
//////
//////            randomDeath();
//////        }catch(MyException me){
//////            me.showMessageBox();
//////        }
//////        catch(Exception e){
//////            throw e;
//////        }
//////    }

//////    private void getNovelty(Chromosome chrm, ArrayList<Chromosome> archive){
//////        ArrayList<Chromosome> entirePopulation = new ArrayList<Chromosome>();
//////        ArrayList<Double> diverse = new ArrayList<Double>();
//////
//////        double distSqrMean = 0;
//////        double roMin;
//////
//////        entirePopulation.addAll(archive);
//////        entirePopulation.addAll(chromosomes_);
//////
//////        double[] mean = new double[entirePopulation.size()];
//////
//////        if (archive.size()<ARCHIVE_MAX){
//////            //get mean
//////            for (int ofsp = 0; ofsp < userInput_.totalDecisionVars; ofsp++) {
//////                diverse.clear();
//////                for (int k = 0; k < entirePopulation.size(); k++) {
//////                    diverse.add(entirePopulation.get(k).vals.get(ofsp));
//////                }
//////                mean[ofsp] = MyMath.mean(diverse);
//////            }
//////
//////            //now get variance - actually ofsp amusing average mean square distane
//////            for (int ofsp = 0; ofsp < entirePopulation.size(); ofsp++) {
//////                distSqrMean += Math.pow(MyMath.norm(entirePopulation.get(ofsp).vals, diverse),2);
//////            }
//////            distSqrMean = distSqrMean/entirePopulation.size();
////////            roMin  = distSqrMean;
//////            roMin = 0;
//////        }
//////        else{
//////
//////        }
//////    }


    
    
    
    
    
//////    /**
//////     * Seems working 
//////     * 
//////     */
//////    private void duplicateSatisfactionSortAndReplace(){
//////        Chromosome tempChrome;
//////        ArrayList<Double> tempVals;
//////        int newSize;
//////        int tempSize;
//////        ArrayList<Integer> randIdx;
//////
//////        for (Chromosome chrm : chromosomes_) {
//////            chrm.sortBy = userInput_.solutionBy; //MUST DO before sorting
//////        }
////////        for (Chromosome chrm : chromosomes_) {
////////            chrm.sortBy = Chromosome.BY_RO; //MUST DO before sorting
////////            try{
////////            chrm.tempRo = getRoValue(chrm);
////////            }catch (Exception e){
////////                e.printStackTrace();
////////            }
////////        }
//////        
//////        bStagnant = false;
//////        if(isStagnant()){
//////            Collections.sort(chromosomes_);
//////            bStagnant = true;
//////            System.diverse.println("*** diverse -- removed...");
//////            tempSize = (int)(3*SAME_BEST_VAL_PERCENT*userInput_.population);
//////
//////            
//////            for (int ofsp= 0; ofsp < tempSize; ofsp++) {
//////                //chromosomes_.remove(0);        
//////                tempVals = chromosomes_.get(ofsp).negateVals();                
//////                chromosomes_.get(ofsp).setVals(tempVals);
//////            }
//////            
////////            try{
////////                //System.diverse.println("be " + chromosomes_.get(0).getValsCopy());
////////                mutation(new ArrayList<Chromosome>(chromosomes_.subList(0, tempSize)));
////////            }catch(UnsupportedDataTypeException udte){
////////                System.diverse.println(udte.getLocalizedMessage());
////////                Application.getInstance().exit();
////////            }
//////        }
//////
//////       
//////        ArrayList<Chromosome> newRandPop = new ArrayList<Chromosome>();
//////
//////
//////
//////         
//////        if(chromosomes_.size() < userInput_.population){            
//////            newSize = userInput_.population-chromosomes_.size();
//////            //randomly make copy of existing chromosomes...
//////            //randIdx = MyRandom.randperm(0, chromosomes_.size()-1);
//////
//////            try {
//////                initializeChromosomes(newRandPop, newSize);
//////            } catch (Exception e) {
//////                e.printStackTrace();
//////                Application.getInstance().exit();
//////            }
//////            
//////            for (int ofsp = 0; ofsp < newSize; ofsp++) {
//////                chromosomes_.addAll(newRandPop);
//////            }
//////
//////            Collections.sort(chromosomes_);
//////            
//////        }else{            
//////            Collections.sort(chromosomes_);
//////            
//////            for (int ofsp = chromosomes_.size()-1; ofsp >= 0; ofsp--) {
//////                if(chromosomes_.size() == userInput_.population)
//////                    break;
//////
//////                tempChrome = chromosomes_.get(ofsp);            
//////                for (Chromosome chrm : chromosomes_) {
//////                    //NOTE THIS STEP .... you can use the commented one as well.....
//////                    //if(chrm != childChrome && chrm.getRankComponents().containsAll(childChrome.getRankComponents())){
//////                    if(chrm != tempChrome && chrm.getValsCopy().containsAll(tempChrome.getValsCopy())){
//////                        chromosomes_.remove(ofsp);                    
//////                        break;
//////                    }
//////                }
//////            }
//////
//////            if(chromosomes_.size() > userInput_.population){ //still more...
//////                chromosomes_ = new ArrayList<Chromosome>(chromosomes_.subList(0, userInput_.population));//pick best ones.
//////            }
//////        }
////////        int d;
////////        ArrayList<Chromosome> newRandPop = new ArrayList<Chromosome>();
////////        d = (int)Math.round(this.REPLACE_PERCENT*userInput_.population);
////////        
////////        try {
////////            initializeChromosomes(newRandPop, d);  
////////        } catch (Exception e) {
////////            e.printStackTrace();
////////            Application.getInstance().exit();
////////        }
////////
////////        for (int ofsp = 0; ofsp < d; ofsp++) {
////////            chromosomes_.set(userInput_.population-1-ofsp, newRandPop.get(ofsp));
////////        } 
//////            randomDeath();
//////    }
    
// </editor-fold>    
     
} //End of class definition
