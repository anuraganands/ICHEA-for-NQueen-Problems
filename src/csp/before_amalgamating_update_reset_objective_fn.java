///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//
//package csp;
//
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.LinkedHashSet;
//import java.util.Set;
//import org.jdesktop.application.Application;
//
//
//
///**
// *
// * @author Anurag
// */
//public class Chromosome implements Comparable, Cloneable {
//    public static final int BY_FITNESS = 1;
//    public static final int BY_VIOLATIONS = 2;
//    public static final int BY_SATISFACTIONS = 3;
//    public static final int BY_RO = 4;
//    public static final int BY_IMMUNITY = 5;
//        
//    private ArrayList<Double> vals_;
//    private ArrayList<Integer> position_;
//    private ArrayList<Double> constFx_;//constraint functions
//    private ArrayList<Double> fitness_;
//    private ArrayList<Integer> violations_;
//    private ArrayList<Integer> satisfactions_;
//    private UserInput userInput_;
//    private ExternalData externalData_;
//    private int noProgressCounter;
//    private int immunity;
//
//    public boolean isStagnant(final int noProgressLimit) {
//        boolean isStagnant;
//        if(noProgressCounter>noProgressLimit)
//            isStagnant = true;
//        else
//            isStagnant = false;
//        
//        return isStagnant;
//    }
//
////    public int fitnessRank;
//    public Double tempRo;
//    public int sortBy;
//    public static int totalEvals=0;
//    //private int rankingType;
//
//    private Chromosome(){
//        vals_ = new ArrayList<Double>();
//        position_ = new ArrayList<Integer>();
//        constFx_ = new ArrayList<Double>();
//        fitness_ = new ArrayList<Double>();
//        violations_ = new ArrayList<Integer>();
//        satisfactions_ = new ArrayList<Integer>();
//        tempRo = -1.0; //Invalid negative value
//        sortBy = BY_VIOLATIONS; //default sort option
//        noProgressCounter = 0;
//        immunity = 0;
//        //rankingType = BY_VIOLATIONS;
//    }
//
//    public Chromosome(UserInput userInput){
//        this();
//        userInput_ = userInput;
//        externalData_ = null;
//        if(userInput_ == null){
//            System.err.println("No user input provided.");
//            Application.getInstance().exit();
//        }
////        fitnessRank = Integer.MAX_VALUE;
//    }
//    public Chromosome(int sortValue, UserInput userInput){
//        this(userInput);
//        sortBy = sortValue; //default sort option
//    }
//
//    public Chromosome(ExternalData externalData){
//        this();
//        externalData_ = externalData;
//        userInput_ = externalData_.getUserInput();
//
//        if(userInput_ == null || this.externalData_ == null){
//            System.err.println("No user input provided or empty external data.");
//            Application.getInstance().exit();
//        }
//    }
//    
//    public Chromosome(int sortValue, ExternalData externalData){
//        this(externalData);
//        sortBy = sortValue; //default sort option
//    }
//
//    public ArrayList<Double> negateVals(){
//        ArrayList<Double> negVals = new ArrayList<Double>();  
//        double val;
//        for (int i = 0; i < vals_.size(); i++) {
//            val = userInput_.minVals.get(i)+userInput_.maxVals.get(i) - vals_.get(i);
//            if(val<userInput_.minVals.get(i))
//                val = userInput_.minVals.get(i);
//            if(val>userInput_.maxVals.get(i))
//                val = userInput_.maxVals.get(i);
//            
//            negVals.add(val);
//        }
//        return negVals;
//    }
//    
//    /**
//     * The integer value denoting the rank of an individual in a population.
//     * the lower the rank the better the value. It has range from [0 Total_Constraints-1]
//     * Gives total violations or total satisfaction in terms of rank<br>
//     * for violation <br>
//     * rank = size of violation<br>
//     * for satisfaction <br>
//     * rank = total constraint - size of satisfaction <br>
//     * Other types of rankings are NOT implemented
//     * @return returns rank - values from good to worse. lower value indicate better chromosome.
//     * @throws UnsupportedOperationException
//     */
//    public int getRank() throws UnsupportedOperationException{
//        int rank;
//        if(sortBy == BY_VIOLATIONS){
//            rank = violations_.size(); //the higher the violation the lower the rank (ascending order from 0...maxval)
//        }else if(sortBy == BY_SATISFACTIONS){
//            rank = userInput_.totalConstraints - satisfactions_.size();
//        }else{
//            rank = -1;
//            throw new UnsupportedOperationException("Incorrect Ranking Type Defined.");
//        }
//        return rank;
//    }
//
//    /**
//     * Depending on the type of rank used in calculation it return the numeric
//     * values of the rank components
//     * <BR>
//     * In case of violation - it return all violated constraint set
//     * <BR>
//     * In case of satisfaction - it return all satisfied constraint set
//     * @return 
//     */
//    public ArrayList<Integer> getRankComponents(){
//        ArrayList<Integer> rankComponents;
//
//        if(sortBy == BY_VIOLATIONS){
//            rankComponents = violations_; //the higher the violation the lower the rank (ascending order from 0...maxval)
//        }else if(sortBy == BY_SATISFACTIONS){
//            rankComponents = satisfactions_;
//        }else{
//            rankComponents = null;
//            throw new UnsupportedOperationException("Incorrect Ranking Type Defined.");
//        }
//        return rankComponents;
//    }
//    
//    /**
//     * Checks if the child is acceptable in inter-marriage crossover.
//     * It is possible that child may not contain any trait of a parent. In this
//     * case it is necessary to check for genuineness with this function.
//     * If constraint violation is to be checked then <BR>
//     * child must not violate any other constraint than its parent's already 
//     * violated constraint. <BR>
//     * If constraint satisfaction is to be checked then <BR>
//     * child must satisfy same or more constraints.
//     * @param child
//     * @return 
//     */
//    public boolean isMyChild(Chromosome child){
//        boolean result;
//        
//        if(sortBy == BY_VIOLATIONS){ //child must not violate any other constraint than its parent's already violated constraint.
//            result = this.violations_.containsAll(child.violations_);
//            if(result){
//                for (int i = 0; i < position_.size(); i++) {
//                    if(position_.get(i) != child.position_.get(i)){
//                        result = false;
//                        break;
//                    }
//                }
//                    
//            }
//        }else if(sortBy == BY_SATISFACTIONS){ //child must satisfy same or more constraints.
//            result = child.satisfactions_.containsAll(this.satisfactions_);
//        }else{
//            throw new UnsupportedOperationException("Incorrect Ranking Type Defined.");
//        }   
//        
//        return result;
//    }
//    
//    
//    /**
//     * It remove duplicate values - It is problem dependant.
//     */
//    private void repairVal(){
//        Set<Double> s = new LinkedHashSet<Double>(this.vals_);
//        this.vals_ = new ArrayList<Double>(s);      
//    }
//    
//    /**
//     * This is tricky one if we deal with both numeric and ordinal data<BR>
//     * - for numeric data vals are features/dimensions where the total dimension 
//     * remain constant. as in my circle example. val size will remain 2.<BR>
//     * - for ordinal data vals are just constraint representation. it can keep
//     * on growing for any size (in some problems up to max satisfaction)
//     * @param vals - 
//     * <BR> - for numeric data - val is dimension/feature value.
//     * <BR> - for ordinal data - val is constraint representation.
//     */
//    public void appendVal(Double vals) {
//        this.vals_.add(vals);
//
//        if(this.vals_.size() == this.userInput_.totalDecisionVars){
//            updateObjectiveFunctionVars();
//        }
//    }
//
//    /**
//     * Updates the vals with the new vals_ from argument
//     * @param vals_ to be assigned to chromosome
//     */
//    public void setVals(ArrayList<Double> vals_) {
//        this.vals_ = vals_;
//        resetObjectiveFunctionVars();
//    }
//
//    public void replaceVal(int idx, Double val){
//        this.vals_.set(idx, val);
//        resetObjectiveFunctionVars();
//    }
//
//    /**
//     * returns the value specified by the index
//     * @param index the index of the decision value
//     * @return the value of specified index of the decision value array
//     */
//    public Double getVals(int index){
//        return this.vals_.get(index);
//    }
//
//    /**
//     * returns the value array;
//     * @return the value array
//     */
//    public ArrayList<Double> getVals() {
//        return this.vals_;
//    }
//
//    /**
//     * Checks to see if there is any violation
//     * if rank is 0, it has no violation.
//     * @return true or false for violation
//     */
//    public boolean isSolution(){
//        return (getRank() == 0);
//    }
//
//    /**
//     * Checks to see if the input chromosome has the same violation as the
//     * current object.
//     * @param from the chomosome with whom comparision is to be made.
//     * @return true or false for same violation
//     */
//    public boolean hasSameRankComponent(Chromosome from){
//        boolean bsame;        
//        if(sortBy == BY_VIOLATIONS){
//            if(this.violations_.size() == this.userInput_.totalConstraints)
//                bsame = false;
//            else            
//                bsame = this.violations_.containsAll(from.violations_) || from.violations_.containsAll(this.violations_); //the higher the violation the lower the rank (ascending order from 0...maxval)
//        }else if(sortBy == BY_SATISFACTIONS){
//            if(this.satisfactions_.isEmpty() || from.satisfactions_.isEmpty())
//                bsame = false;
//            else
//                bsame = this.satisfactions_.containsAll(from.satisfactions_) || from.satisfactions_.containsAll(this.satisfactions_);
//        }else{
//            throw new UnsupportedOperationException("Incorrect Ranking Type Defined.");
//        }
//        return bsame;
//    }
//
//    public ArrayList<Integer> getPositionRelConst() {
//        return position_;
//    }
//  
//    
//
//    
//
//    /**
//     * NOT TESTED... TEST IT FIRST BEFORE USING IT
//     */
//    private void updateObjectiveFunctionVars(){
//        //DO NOT DELETE THESE LINES
//        //<<
//            totalEvals++;
//            int prevRank;
//            int curRank;
//            prevRank = this.getRank();
//        //>>
//
//        //TODO... Call your objective function here....
//        //Start<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//            objFunction_NQueenII();
//        //End>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//                
//        //DO NOT DELETE THESE LINES
//        //<<
//            curRank = this.getRank();        
//            if(curRank>=prevRank){ //worse or same => stagnant
//                noProgressCounter++;
//            }else{
//                noProgressCounter = 0;
//            }
//        //>>
//    }
//    
//    /**
//     * NOT TESTED.... TEST IT PROPERLY BEFORE USE
//     */
//    private void resetObjectiveFunctionVars(){
//        totalEvals++;
//        //DO NOT DELETE THESE LINES
//        //<<
//            int prevRank;
//            int curRank;
//            prevRank = this.getRank();
//        //>>
//
//       
//        //TODO... Call your objective function here....
//        //Start<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//            objFunction_NQueenII();
//        //End>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//            
//        //DO NOT DELETE THESE LINES
//        //<<
//            curRank = this.getRank();        
//            if(curRank>=prevRank){ //worse or same => stagnant
//                noProgressCounter++;
//            }else{
//                noProgressCounter = 0;
//            }
//        //>>
//    }
//
//    /**
//     * checks if marriage is compatible in a sense that two chromosomes should
//     * have at least two opposite positions which can help in finding the point
//     * of intersection. for example {-1,-1, 0 1} and {-1, 1, 0 , -1} has opposite
//     * positions at 2nd and 4th place.
//     * 
//     * @param c another couple
//     * @return 
//     */
//    public boolean  isMarriageCompatible(Chromosome c){        
//        int count = 0;
//        boolean result = false;
//        
//        if(this.sortBy != BY_VIOLATIONS){
//            System.err.println("this.sortBy must be BY_VIOLATION. \n Error in isMarriageCompatible()");
//            Application.getInstance().exit();
//        }
//        
//        CspProcess.bringCloserRatio = 0.5;
//        
//        if(this.hasSameRankComponent(c) && this.getRank() != userInput_.totalConstraints-1){
//            CspProcess.bringCloserRatio = 0.10;
//            result = true;            
//        }else{ 
//            CspProcess.bringCloserRatio = 0.5;
//            //for (int v = 1; v <= userInput_.totalConstraints; v++) {
//            for (Integer v : violations_) {
//                if(this.position_.get(v-1) != c.position_.get(v-1)){
//                    count++;
//                }
//                
////                CspProcess.abToHoJaFlag = false;
////                if((violations_.size() == 1 && count == 1)){
////                   immunity+=100;
////                   CspProcess.abToHoJaFlag = true;
////                   System.out.println("ab to ho ja.");
////                }
////                
////                if((c.violations_.size() == 1 && count == 1)){
////                   c.immunity+=100;
////                   CspProcess.abToHoJaFlag = true;                   
////                }
//                
////                if(count == 2) {// || (violations_.size() == 1 && count == 1) || (c.violations_.size() == 1 && count == 1)){
////                    immunity+=2;
////                    c.immunity+=2;
////                    result = true;
////                    break;
////                }
//            }
//            
//                if(count > 0) {// || (violations_.size() == 1 && count == 1) || (c.violations_.size() == 1 && count == 1)){
//                    immunity+=count;
//                    c.immunity+=count;
//                    result = true;
//                }
//        }
//        
//        return  result;
//    }
//
//    public void useImmunity() {
//        this.immunity--;
//        if(immunity < 0){
//            immunity = 0;
//        }
//    }
//
//    public int getImmunity() {
//        return immunity;
//    }
//    
//    
//    
//    
//    private void objFunction_deb(){
//    //<< deb problem
//
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        this.fitness_.add(this.vals_.get(0));
//        this.fitness_.add((1+this.vals_.get(1))/this.vals_.get(0));
//
//        ArrayList<Double> x = this.vals_;
//        ArrayList<Double> f = this.fitness_;
//
//        this.violations_.clear();
//        this.constFx_.clear();
//        this.position_.clear();
//        
//        //should bes...
//        this.constFx_.add(x.get(1) + 9*x.get(0)-6.0); //=0
////        this.constFx_.add(x.get(1) + 9*x.get(0)-6.1); //=0
//        this.constFx_.add(-x.get(1) + 9*x.get(0)-1.0); //=0
////        this.constFx_.add(-x.get(1) + 9*x.get(0)-1.1); //=0
//        this.constFx_.add(x.get(0)-1); //=0
//        this.constFx_.add(x.get(0)); //=0
//        
//        if(constFx_.size() != this.userInput_.totalConstraints){
//            System.err.println("constFx size is not same as total constraints");
//            Application.getInstance().exit();
//        }
//        
//        for (Double cf : this.constFx_) {
//            if(cf<0){
//                position_.add(-1);
//            }else if (cf>0){
//                position_.add(1);
//            }else{
//                position_.add(0);
//            }            
//        }
//        
//        //if(x.get(1) + 5.3*x.get(0) < 6){
//        if(this.constFx_.get(0) < 0){
//            this.violations_.add(1);
//        }
//
//        //if(-x.get(1) + 1.8*x.get(0) < 1){
//        if(this.constFx_.get(1) < 0){
//            this.violations_.add(2);
//        }
//
//        if(this.constFx_.get(2) > 0){
//            this.violations_.add(3);
//        }
//        
//        if(this.constFx_.get(3) < 0){
//            this.violations_.add(4);
//        }
//
////        if(this.constFx_.get(0)<0){
////            this.violations_.add(1);
////        }
////
////        if(this.constFx_.get(1)>0){
////            this.violations_.add(2);
////        }
////
////        if(this.constFx_.get(2)<0){
////            this.violations_.add(3);
////        }
////        
////        if(this.constFx_.get(3)>0){
////            this.violations_.add(4);
////        }
//    //>> deb problem
//    }
//    
//    
//     private void objFunction_imod(){
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        this.fitness_.add(this.vals_.get(0));
//        this.fitness_.add((1+this.vals_.get(1))/this.vals_.get(0));
//
//        ArrayList<Double> x = this.vals_;
//        ArrayList<Double> f = this.fitness_;
//
//        this.violations_.clear();
//        this.constFx_.clear();
//        this.position_.clear();  
//        double dt = 0.1;
//        
//        //vars = 8
//        //const = 8
//        
//       
//        
//        this.constFx_.add(x.get(0) - 0.25428722 - 0.18324757*x.get(3)*x.get(2)*x.get(8));
//        this.constFx_.add(x.get(1) - 0.37842197 - 0.16275449*x.get(0)*x.get(9)*x.get(5));
//        this.constFx_.add(x.get(2) - 0.27162577 - 0.16955071*x.get(0)*x.get(1)*x.get(9));
//        this.constFx_.add(x.get(3) - 0.19807914 - 0.15585316*x.get(6)*x.get(0)*x.get(5));
//        this.constFx_.add(x.get(4) - 0.44166728 - 0.19950920*x.get(6)*x.get(5)*x.get(2));
//        this.constFx_.add(x.get(5) - 0.14654113 - 0.18922793*x.get(7)*x.get(4)*x.get(9));
//        this.constFx_.add(x.get(6) - 0.42937161 - 0.21180484*x.get(1)*x.get(4)*x.get(7));
//        this.constFx_.add(x.get(7) - 0.07056438 - 0.17081208*x.get(0)*x.get(6)*x.get(5));
//        this.constFx_.add(x.get(8) - 0.34504906 - 0.19612740*x.get(9)*x.get(5)*x.get(7));
//        this.constFx_.add(x.get(9) - 0.42651102 - 0.21466544*x.get(3)*x.get(7)*x.get(0));
//
//
//
//
//        
//        if(constFx_.size() != this.userInput_.totalConstraints){
//            System.err.println("constFx size is not same as total constraints");
//            Application.getInstance().exit();
//        }
//        
//        for (Double cf : this.constFx_) {
//            if(cf<0){
//                position_.add(-1);
//            }else if (cf>0){
//                position_.add(1);
//            }else{
//                position_.add(0);
//            }            
//        }
//        
//        
//        if(this.constFx_.get(0)>dt || this.constFx_.get(0)< -dt){
//            this.violations_.add(1);
//        }
//
//        if(this.constFx_.get(1)>dt|| this.constFx_.get(1)< -dt){
//            this.violations_.add(2);
//        }
//
//        if(this.constFx_.get(2)> dt || this.constFx_.get(2)<-dt){
//            this.violations_.add(3);
//        }
//        
//        if(this.constFx_.get(3)>dt || this.constFx_.get(3)< -dt){
//            this.violations_.add(4);
//        }
//
//        if(this.constFx_.get(4)>dt|| this.constFx_.get(4)< -dt){
//            this.violations_.add(5);
//        }
//
//        if(this.constFx_.get(5)> dt || this.constFx_.get(5)<-dt){
//            this.violations_.add(6);
//        }
//                
//        if(this.constFx_.get(6)>dt || this.constFx_.get(6)< -dt){
//            this.violations_.add(7);
//        }  
//        
//        if(this.constFx_.get(7)>dt || this.constFx_.get(7)< -dt){
//            this.violations_.add(8);
//        }  
//         
//        if(this.constFx_.get(8)>dt || this.constFx_.get(8)< -dt){
//            this.violations_.add(9);
//        }  
//          
//        if(this.constFx_.get(9)>dt || this.constFx_.get(9)< -dt){
//            this.violations_.add(10);
//        }             
//     }
//    
//    //<<My CIRCLEs Test Problem
//    private void objFunction_circles(){
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        this.fitness_.add(this.vals_.get(0));
//        this.fitness_.add((1+this.vals_.get(1))/this.vals_.get(0));
//
//        ArrayList<Double> x = this.vals_;
//        ArrayList<Double> f = this.fitness_;
//
//        this.violations_.clear();
//        this.constFx_.clear();
//        this.position_.clear();       
//        
//        //should bes...
//        double r1 = 2.0;
//        double r2 = 1.0;
//        double r3 = 1.0;
//        double outerBoundary = 10.0;
//        
//        this.constFx_.add((x.get(1)-2)*(x.get(1)-2) + (x.get(0)-2)*(x.get(0)-2)-r1*r1); //=0
//        this.constFx_.add((x.get(1)-2)*(x.get(1)-2) + Math.pow(x.get(0)-4.90, 2)-r2*r2); //=0
//        this.constFx_.add((x.get(1)-1)*(x.get(1)-1) + (x.get(0)-4.50)*(x.get(0)-4.50)-r3*r3); //=0
//        this.constFx_.add(x.get(1)-outerBoundary);
//        this.constFx_.add(x.get(1)+outerBoundary);
//        this.constFx_.add(x.get(0)-outerBoundary);
//        this.constFx_.add(x.get(0)+outerBoundary);
//        //this.constFx_.add(x.get(1) - x.get(0));
//        
//        if(constFx_.size() != this.userInput_.totalConstraints){
//            System.err.println("constFx size is not same as total constraints");
//            Application.getInstance().exit();
//        }
//        
//        for (Double cf : this.constFx_) {
//            if(cf<0){
//                position_.add(-1);
//            }else if (cf>0){
//                position_.add(1);
//            }else{
//                position_.add(0);
//            }            
//        }
//   
//        //for circles
//        if(this.constFx_.get(0)>0){
//            this.violations_.add(1);
//        }
//
//        if(this.constFx_.get(1)>0){
//            this.violations_.add(2);
//        }
//        
//        if(this.constFx_.get(2)>0){
//            this.violations_.add(3);
//        } 
//        
//        
//        //for outer boundaries.
//        if(this.constFx_.get(3)>0){
//            this.violations_.add(4);
//        } 
//        if(this.constFx_.get(4)<0){
//            this.violations_.add(5);
//        }
//        if(this.constFx_.get(5)>0){
//            this.violations_.add(6);
//        }
//        if(this.constFx_.get(6)<0){
//            this.violations_.add(7);
//        } 
//        
//        
////////        if(this.constFx_.get(7)>0){
////////            this.violations_.add(8);
////////        } 
//    }
//    
//    
//    
//    
//    
//    
//    private void objFunction_circumscribedCircle(){
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        this.fitness_.add(this.vals_.get(0));
//        this.fitness_.add((1+this.vals_.get(1))/this.vals_.get(0));
//
//        ArrayList<Double> x = this.vals_;
//        ArrayList<Double> f = this.fitness_;
//
//        this.violations_.clear();
//        this.constFx_.clear();
//        this.position_.clear();
//        
//        //should bes...
//        double br = 4.0;
//        double sr = 0.1;
//        
//        this.constFx_.add(x.get(1)*x.get(1) + x.get(0)*x.get(0)-br*br); //=0
//        this.constFx_.add(x.get(1)*x.get(1) + x.get(0)*x.get(0)-sr*sr); //=0
//        
//        if(constFx_.size() != this.userInput_.totalConstraints){
//            System.err.println("constFx size is not same as total constraints");
//            Application.getInstance().exit();
//        }
//        
//        for (Double cf : this.constFx_) {
//            if(cf<0){
//                position_.add(-1);
//            }else if (cf>0){
//                position_.add(1);
//            }else{
//                position_.add(0);
//            }            
//        }
//   
//        if(this.constFx_.get(0)>0){
//            this.violations_.add(1);
//        }
//
//        if(this.constFx_.get(1)>0){
//            this.violations_.add(2);
//        }
//       
//    //>> deb problem
//    }
//    
//    
////    private void objFunction_test(){
////        
////        this.fitness_ = new ArrayList<Double>();
////        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
////
////        this.fitness_.add(this.vals_.get(0));
////        this.fitness_.add((1+this.vals_.get(1))/this.vals_.get(0));
////
////        ArrayList<Double> x = this.vals_;
////        ArrayList<Double> f = this.fitness_;
////
////        this.violations_.clear();
////        this.constFx_.clear();
////        this.position_.clear();
////        
////        double best_val_found = 0.2415051288;
////	double eps = 1;
////	double dt = 0.001;
////        double s2 = Math.sqrt(2.0);
////        //var
//////        x {1..5} >= 0, <= 4;
////                //cons = 5
////		//var = 5
////        this.constFx_.add( best_val_found + eps - (2*x.get(0)*(x.get(0)-x.get(1)-1) + 1 + Math.pow(x.get(1),2) + Math.pow(x.get(2) - 1, 2) + Math.pow(x.get(3) - 1.0,4.0) + Math.pow(x.get(4) - 1.0,6.0)));
////        this.constFx_.add(dt - (Math.pow(x.get(0),2) * x.get(3) + Math.sin(x.get(3) - x.get(4)) - 2.0 * s2)); //=0
////        this.constFx_.add(dt + Math.pow(x.get(0),2) * x.get(3) + Math.sin(x.get(3) - x.get(4)) - 2.0 * s2); //=0
////        this.constFx_.add(dt - (x.get(1) + Math.pow(x.get(2)*1.0,4.0)*Math.pow(x.get(3),2) - 8 - s2)); //=0
////        this.constFx_.add(dt +  x.get(1) + Math.pow(x.get(2)*1.0,4.0)*Math.pow(x.get(3),2) - 8 - s2); //=0
////
////        if(constFx_.size() != this.userInput_.totalConstraints){
////            System.err.println("constFx size is not same as total constraints");
////            Application.getInstance().exit();
////        }
////        
////        for (Double cf : this.constFx_) {
////            if(cf<0){
////                position_.add(-1);
////            }else if (cf>0){
////                position_.add(1);
////            }else{
////                position_.add(0);
////            }            
////        }
////
////        
////        if(this.constFx_.get(0)<0){
////            this.violations_.add(1);
////        }
////
////        if(this.constFx_.get(1)<0){
////            this.violations_.add(2);
////        }
////        
////        if(this.constFx_.get(2)<0){
////            this.violations_.add(3);
////        }
////
////        if(this.constFx_.get(3)<0){
////            this.violations_.add(4);
////        }
////        if(this.constFx_.get(4)<0){
////            this.violations_.add(5);
////        }
////       
////        
////    }
//    
//        private void objFunction_test(){
//        
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        this.fitness_.add(this.vals_.get(0));
//        this.fitness_.add((1+this.vals_.get(1))/this.vals_.get(0));
//
//        ArrayList<Double> x = this.vals_;
//        ArrayList<Double> f = this.fitness_;
//
//        this.violations_.clear();
//        this.constFx_.clear();
//        this.position_.clear();
//        
//        double best_val_found = 0.2415051288;
//	double eps = 1;
//	double dt = 0.1;
//        double s2 = Math.sqrt(2.0);
//        //var
////        x {1..5} >= 0, <= 4;
//                //cons = 5
//		//var = 5
//        this.constFx_.add( best_val_found + eps - (2*x.get(0)*(x.get(0)-x.get(1)-1) + 1 + Math.pow(x.get(1),2) + Math.pow(x.get(2) - 1, 2) + Math.pow(x.get(3) - 1.0,4.0) + Math.pow(x.get(4) - 1.0,6.0)));
//        this.constFx_.add(Math.pow(x.get(0),2) * x.get(3) + Math.sin(x.get(3) - x.get(4)) - 2.0 * s2); //=0
//        this.constFx_.add(x.get(1) + Math.pow(x.get(2)*1.0,4.0)*Math.pow(x.get(3),2) - 8 - s2); //=0
//        		
//
//        if(constFx_.size() != this.userInput_.totalConstraints){
//            System.err.println("constFx size is not same as total constraints");
//            Application.getInstance().exit();
//        }
//        
//        for (Double cf : this.constFx_) {
//            if(cf<0){
//                position_.add(-1);
//            }else if (cf>0){
//                position_.add(1);
//            }else{
//                position_.add(0);
//            }            
//        }
//
//        
//        if(this.constFx_.get(0)<0){
//            this.violations_.add(1);
//        }
//
//        if(this.constFx_.get(1)>dt|| this.constFx_.get(1)< -dt){
//            this.violations_.add(2);
//        }
//
//        if(this.constFx_.get(2)> dt || this.constFx_.get(2)<-dt){
//            this.violations_.add(3);
//        }
//       
//        
//    }
//        
//     private void objFunction_bangleCircle(){
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        this.fitness_.add(this.vals_.get(0));
//        this.fitness_.add((1+this.vals_.get(1))/this.vals_.get(0));
//
//        ArrayList<Double> x = this.vals_;
//        ArrayList<Double> f = this.fitness_;
//
//        this.violations_.clear();
//        this.constFx_.clear();
//        this.position_.clear();
//        
//        //should bes...
//        double br = 5.0;
//        double sr = 4.9;
//        double outerBoundary = 30.0;
//        
//        this.constFx_.add(x.get(1)*x.get(1) + x.get(0)*x.get(0)-br*br); //=0
//        this.constFx_.add(x.get(1)*x.get(1) + x.get(0)*x.get(0)-sr*sr); //=0
//        this.constFx_.add(x.get(1)*x.get(1) + Math.pow(x.get(0)+2*br-1.0,2)-br*br); //=0
//        this.constFx_.add(x.get(1)*x.get(1) + Math.pow(x.get(0)+2*br-1.0,2)-sr*sr); //=0
//        
//        this.constFx_.add(x.get(1)-outerBoundary);
//        this.constFx_.add(x.get(1)+outerBoundary);
//        this.constFx_.add(x.get(0)-outerBoundary);
//        this.constFx_.add(x.get(0)+outerBoundary);
//        
//     
//        if(constFx_.size() != this.userInput_.totalConstraints){
//            System.err.println("constFx size is not same as total constraints");
//            Application.getInstance().exit();
//        }
//        
//        for (Double cf : this.constFx_) {
//            if(cf<0){
//                position_.add(-1);
//            }else if (cf>0){
//                position_.add(1);
//            }else{
//                position_.add(0);
//            }            
//        }
//   
//        if(this.constFx_.get(0)>0){
//            this.violations_.add(1);
//        }
//
//        if(this.constFx_.get(1)<0){
//            this.violations_.add(2);
//        }
//        
//        if(this.constFx_.get(2)>0){
//            this.violations_.add(3);
//        }
//
//        if(this.constFx_.get(3)<0){
//            this.violations_.add(4);
//        }
//       
//        //for outer boundaries.
//        if(this.constFx_.get(4)>0){
//            this.violations_.add(5);
//        } 
//        if(this.constFx_.get(5)<0){
//            this.violations_.add(6);
//        }
//        if(this.constFx_.get(6)>0){
//            this.violations_.add(7);
//        }
//        if(this.constFx_.get(7)<0){
//            this.violations_.add(8);
//        } 
//        
//    }
//     
//    
//    private void objFunction_broyden10(){
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//        
//        //no fitness function
//        
//        //ArrayList<Double> x = this.vals_;
//        Double x[] = new Double[this.userInput_.totalConstraints];
//        this.vals_.toArray(x);
//        
//        this.violations_.clear();
//        
//        if(x[0]*(2+5*x[0]*x[0]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1])) != 0){
//            this.violations_.add(1);
//        }
//                
//        if(x[1]*(2+5*x[1]*x[1]) + 1 - (x[0]*(1+x[0]) + x[2]*(1+x[2])) != 0){
//            this.violations_.add(2);
//        }
//        
//        if(x[2]*(2+5*x[2]*x[2]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[3]*(1+x[3])) != 0){
//            this.violations_.add(3);
//        }
//        
//        if(x[3]*(2+5*x[3]*x[3]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[4]*(1+x[4])) != 0){
//            this.violations_.add(4);
//        }
//        
//        if(x[4]*(2+5*x[4]*x[4]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[3]*(1+x[3])+ x[5]*(1+x[5])) != 0){
//            this.violations_.add(5);
//        }
//
//        if(x[5]*(2+5*x[5]*x[5]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[3]*(1+x[3])+ x[4]*(1+x[4]) + x[6]*(1+x[6])) != 0){
//            this.violations_.add(6);
//        }
//
//        if(x[6]*(2+5*x[6]*x[6]) + 1 - (x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[3]*(1+x[3]) + x[4]*(1+x[4])+ x[5]*(1+x[5]) + x[7]*(1+x[7])) != 0){
//            this.violations_.add(7);
//        }
//
//        if(x[7]*(2+5*x[7]*x[7]) + 1 - (x[2]*(1+x[2]) + x[3]*(1+x[3]) + x[4]*(1+x[4]) + x[5]*(1+x[5])+ x[6]*(1+x[6]) + x[8]*(1+x[8])) != 0){
//            this.violations_.add(8);
//        }
//
//        if(x[8]*(2+5*x[8]*x[8]) + 1 - (x[3]*(1+x[3]) + x[4]*(1+x[4]) + x[5]*(1+x[5]) + x[6]*(1+x[6])+ x[7]*(1+x[7]) + x[9]*(1+x[9])) != 0){
//            this.violations_.add(9);
//        }
//
//        if(x[9]*(2+5*x[9]*x[9]) + 1 - (x[4]*(1+x[4]) + x[5]*(1+x[5]) + x[6]*(1+x[6]) + x[7]*(1+x[7])+ x[8]*(1+x[8])) != 0){
//            this.violations_.add(10);
//        }
//        
//    }
//    
//    
//    private void objFunction_biggsc4(){
//    //<< deb problem
//
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        //no fitness function.
//
//        this.violations_.clear();
//        this.constFx_.clear();
//        this.position_.clear();
//        
//        ArrayList<Double> x = this.vals_;
//        double best_val_found = -24.375;
//	double eps = 1; 
//        
//        this.constFx_.add(best_val_found + eps - (-x.get(0)*x.get(2)-x.get(1)*x.get(3)));
//        this.constFx_.add(x.get(0)+x.get(1) -2.5);
//        this.constFx_.add(x.get(0)+x.get(2) -2.5); //=0
//        this.constFx_.add(x.get(0)+x.get(3) -2.5);
//        this.constFx_.add(x.get(1)+x.get(2) -2.0);
//        this.constFx_.add(x.get(1)+x.get(3) -2.0);
//        this.constFx_.add(x.get(2)+x.get(3) -1.5);
//        this.constFx_.add(x.get(0)+x.get(1)+x.get(2)+x.get(3)-5.0);    
//        
//        if(constFx_.size() != this.userInput_.totalConstraints){
//            System.err.println("constFx size is not same as total constraints");
//            Application.getInstance().exit();
//        }
//        
//        for (Double cf : this.constFx_) {
//            if(cf<0){
//                position_.add(-1);
//            }else if (cf>0){
//                position_.add(1);
//            }else{
//                position_.add(0);
//            }            
//        }
//
//        if(this.constFx_.get(0)<0){
//            this.violations_.add(1);
//        }
//        if(this.constFx_.get(1)<0|| this.constFx_.get(1)> 5){
//            this.violations_.add(2);
//        }
//        if(this.constFx_.get(2)<0|| this.constFx_.get(2)> 5){
//            this.violations_.add(3);
//        }
//        if(this.constFx_.get(3)<0|| this.constFx_.get(3)> 5){
//            this.violations_.add(4);
//        }
//        if(this.constFx_.get(4)<0|| this.constFx_.get(4)> 5){
//            this.violations_.add(5);
//        }
//        if(this.constFx_.get(5)<0|| this.constFx_.get(5)> 5){
//            this.violations_.add(6);
//        }
//        if(this.constFx_.get(6)<0|| this.constFx_.get(6)> 5){
//            this.violations_.add(7);
//        }
//        if(this.constFx_.get(7)<0){
//            this.violations_.add(8);
//        }
//    }
//    
//    
//    private void objFunction_h74mod(){
//        
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        //no fitness function.
//
//        this.violations_.clear();
//        Double x[] = new Double[this.userInput_.totalConstraints];
//        this.vals_.toArray(x);
//        
//        double a = 0.55;  
//        double c1 = Math.pow(2.0/3.0,-6);
//
////	cons1 {i in 1..2}: x[i] >= 0;
////cons2 {i in 1..2}: x[i] <= 1200;
////cons3 {i in 3..4}: x[i] >= -a;
////cons4 {i in 3..4}: x[i] <= a;
//
//
//    if(x[0] * (3 + Math.pow(1.0,-6.0)*x[0]*x[0]) + x[1] * (2 + Math.pow(c1*x[1],2.0)) > 5126.49811 + 51.2649811);
//	this.violations_.add(1);
//
//    if(x[3] - x[2] + a <0)
//	this.violations_.add(2);
//
//    if(x[2] - x[3] + a <0)
//	this.violations_.add(3);
//
//    if(1000 * Math.sin(-x[2] - 0.25) + 1000 * Math.sin(-x[3]-0.25) + 894.8 - x[0] != 0)
//	this.violations_.add(4);	
//
//    if(1000 * Math.sin(x[2] - 0.25) + 1000 * Math.sin(x[2]-x[3]-0.25) + 894.8 - x[1] != 0)
//	this.violations_.add(5);
//
//    if(1000 * Math.sin(x[3] - 0.25) + 1000 * Math.sin(x[3]-x[2]-0.25) + 1294.8 != 0)
//	this.violations_.add(6);
//    }
//    
////    public ArrayList<Double> correspondingBoundaryDist(){
////        ArrayList<Double> rank = new ArrayList<Double>();
////        rank.add((double)Math.round(Math.abs(this.vals_.get(1) + 9*this.vals_.get(0) - 6)));
////        rank.add((double)Math.round(Math.abs(-this.vals_.get(1) + 9*this.vals_.get(0) - 1)));
////        rank.add((double)Math.round(Math.abs(-this.vals_.get(0)+1)));
////        
////        return rank;
////    }
//    
//    
//    private void objFunction_StackedGymnasticHalls(){
//        
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        //no fitness function.
//
//        this.violations_.clear();
//        Double x[] = new Double[this.userInput_.totalConstraints];
//        this.vals_.toArray(x);
//        
//        
//        //# Constants
//        double Pi = 3.14159265358979323846;
//        double ZeroPlus = Math.pow(10.0,-12.0);
//
//        //# by dynamics experts
//        double E_steel	= 220000.0;	
//
//        //# by civil engineer
//        double R_f	= 2.0;		
//        double S	= 1.65;	
//
//        //# additional
//        double h_t	= 0.5;		
//        double L_z = 11.5;	
//        ////
//        ////
//        //////#-------------------------------------------------------------------------------
//        //////# 12 variables
//        //////#-------------------------------------------------------------------------------
//        double E = x[0]; //		>=2	,<=10	;#2	# beam spacing [m]
//        double H_s1	= x[1]; //>=5	,<=8	;#7	# height of lower gym [m]
//        double H_s2	=x[2]; //>=5	,<=8	;#7	# height of upper gym [m]
//        double L; //		>=10	,<=50	;#27.4	# span [m]
//        double c = x[3]; //		>=200	,<=500	;#200	# half cover plate width [mm]
//        double d = x[4]; //		>=20	,<=40	;#20	# web thickness [mm]
//        double h_a = x[5]; //		>=800	,<=4000	;#4000	# beam height [mm]
//        double k1; //		>=1	,<=1e8	;#1e7	# modal stiffness [N/m]
//        double m; //		>=1	,<=1e6	;#1000	# modal weight [kg]
//        double P = x[6]; //		>=0	,<=10	;#5	# depth of foundation [m]
//        double q_concrete; //	>=200	,<=2000	;#200	# self weight [kg/m]
//        double t = x[7];// 		>=20	,<=100	;#20	# half cover plate thickness [mm]
//
//
//
//
//        //#ArchitectC01:	
//        L = 27.4;
//        //DynamicsC00:	
//        q_concrete = 200 * Math.pow(E,0.6);
//        //DynamicsC01:	
//        m = q_concrete * E * L / 2.0;
//        //DynamicsC02:	
//        k1 = E_steel/Math.pow(L,3.0)/Math.pow(10.0,6.0)*(3.37*Math.pow(h_a,3.0)*d+42.84*Math.pow(h_a,2.0)*t*c);
//
//
//        //# Constraints by client
//        //ClientC00:	
//        if(H_s1 < 5.5 + ZeroPlus)	//# needed for certain sports
//            this.violations_.add(1);
//        //ClientC01:
//        if(H_s2 < 5.5 + ZeroPlus)	//# needed for certain sports
//            this.violations_.add(2);
//        //ClientC02:	
//        if(E < 2.5 + ZeroPlus)	//# needed to mount sports facilities
//            this.violations_.add(3);
//
//        //# Constraints by contractor
//        //ContractorC00:
//        if(4.5*t > c)
//            this.violations_.add(4);
//        //ContractorC01:
//        if(c > 8.5*t)
//            this.violations_.add(5);
//        //ContractorC02:
//        if(40*d > h_a)
//            this.violations_.add(6);
//        //ContractorC03:
//        if(h_a > 100*d)
//            this.violations_.add(7);
//
//        //# Constraints by dynamics expert
//
//        //DynamicsC03:	
//        if(Math.pow((k1/m),0.5) < 8 * 2 * Pi)
//            this.violations_.add(8);
//
//        //# Constraints by civil engineer
//        //CivilC00:
//        if(0.5 * h_a * d * 235 / Math.pow(3,0.5) < 1000 * L * Math.pow(E,8.0/5)*R_f*S)
//            this.violations_.add(9);
//        //CivilC01:	
//        if(235*t*c*h_a/500.0 < 250 * L * Math.pow(E,8.0/5)*R_f*S)
//            this.violations_.add(10);
//
//        //# Constraints by Geologist
//        //GeologistC00:	
//        if(P > 5 - ZeroPlus)
//            this.violations_.add(11);
//
//        //# Constraints by architect
//        //ArchitectC00:
//        if(H_s1 + H_s2 + h_a/1000.0 + h_t > L_z + P - ZeroPlus)
//            this.violations_.add(12);
//    }
//    
//
//     private void objFunction_TrussDesign01(){
//        
//        this.fitness_ = new ArrayList<Double>();
//        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//
//        //no fitness function.
//
//        this.violations_.clear();
//        Double x[] = new Double[this.userInput_.totalConstraints];
//        this.vals_.toArray(x);
//        
//        //Constants-------------------------------------------------------
//        double E = 210 * Math.pow(10.0, 6.0); //# (Unit = kN/m2, Youngs modulus of steel)
//        double T = 235000.0; //# (Unit = kN/m2, yield stress of steel)
//        double A = 0.25; //# Area of cross section of truss members
//        double r = 0.5; //# Radius of gyration of the cross section of truss members
//        double Pi = 4*Math.atan(1.0);
//        double ZeroPlus = 0.00000001; //# sufficiently small positive number
//
//        double P = 400.0; //#var H >=5, <=10;
//        double H = 6.0; //#var H >=5, <=10;height of truss ???
//        double L = 10.0;
//        
//        double x1, y1;
//        x1 = x[0];
//        y1 = x[1];
//
//        if((P/Math.tan(Math.atan((H-y1)/(L-x1)))) > T*A - ZeroPlus)
//            this.violations_.add(1);
//
//        if(((P*L/H)*Math.tan(Math.atan(y1/x1))) > T*A - ZeroPlus)
//            this.violations_.add(2);
//
//        if((P/Math.sin(Math.atan((H-y1)/(L-x1)))) > (Math.pow(Pi,2)*E/((Math.pow(Math.pow((L-x1)*(L-x1)+(H-y1)*(H-y1),0.5)/r,2)*A - ZeroPlus))))
//            this.violations_.add(3);
//
//
//        if((P/Math.sin(Math.atan((H-y1)/(L-x1)))) > T*A - ZeroPlus)
//            this.violations_.add(4);
//
//        if((P*L/H)/Math.cos(Math.atan(y1/x1)) > Pi*Pi*E/(Math.pow(   Math.pow(x1*x1+y1*y1,0.5)/r,  2)*A - ZeroPlus))
//            this.violations_.add(5);
//
//        if(((P*L/H)/Math.cos(Math.atan(y1/x1))) > T*A - ZeroPlus)
//            this.violations_.add(6);
//
//        if(Math.abs((P*L/H - P/Math.tan(Math.atan((H-y1)/(L-x1))))/Math.cos(Math.atan((H-y1)/x1))) > T*A - ZeroPlus)
//            this.violations_.add(7);
//
//        if ((P*L/H - P/Math.tan(Math.atan((H-y1)/(L-x1))))/Math.cos(Math.atan((H-y1)/x1)) <= 0){
//            if(-((P*L/H - P/Math.tan(Math.atan((H-y1)/(L-x1))))/Math.cos(Math.atan((H-y1)/x1))) > Pi*Pi*E/(Math.pow(   Math.pow(  (x1*x1+(H-y1)*(H-y1)),  0.5)  /r,   2)*A - ZeroPlus))
//                this.violations_.add(8);
//        }
//
//        if(x1 > L - ZeroPlus)
//            this.violations_.add(9);
//
//        if(y1 > H - ZeroPlus)
//            this.violations_.add(10);
//        
//     }
//     
//    private void objFunction_NQueenII(){          
//        this.fitness_.clear();
//        this.fitness_.add(Double.NaN);
//        this.sortBy = Chromosome.BY_SATISFACTIONS;
//        this.satisfactions_.clear();
//
//        repairVal();
//
//        for (int i = 0; i < vals_.size(); i++) {
//            for (int j = i+1; j < vals_.size(); j++) {
//                 if(i!=j){
//                    if(externalData_.isViolated(vals_.get(i).intValue(), vals_.get(j).intValue(),Math.abs(i-j))){
//                        vals_.remove(j);
//                        j--;
//                    }else{
//                        ;
//                    }
//                }
//            }
//        }
//
//        for (Double v : vals_) {
//            satisfactions_.add(v.intValue());
//        }                    
//    }
//    
//    @Override
//    public String toString() {
//        String sp = "  ";
//        String print = "\n"+sp+"{";
//        print += "\n"+sp+sp+"vals:" + vals_.toString();
//        print += "\n"+sp+sp+"Fitness:" + fitness_.toString();
//        print += "\n"+sp+sp+"Violations:" + violations_.toString();
//        print += "\n"+sp+sp+"Satisfaction:" + satisfactions_.toString();
//        print += "\n"+sp+sp+"Rank:" + getRank()+"/[0:"+(userInput_.totalConstraints-1)+"]";
//        print += "\n"+sp+sp+"tempRo:" + tempRo.toString();
//        print += "\n"+sp+"}\n";
//        return print;
//    }
//
//    @Override
//    public int compareTo(Object obj){
//        if (!(obj instanceof Chromosome)) {
//            throw new ClassCastException("Not a Chromosome");
//        }
//        Chromosome c = (Chromosome) obj;
//
//        if(this.sortBy == Chromosome.BY_VIOLATIONS || this.sortBy == Chromosome.BY_SATISFACTIONS) //good to bad - less to more
//            //return this.violations_.size() - c.violations_.size();
//            return this.getRank() - c.getRank();
//        else if(this.sortBy == Chromosome.BY_IMMUNITY){
//            return this.immunity - c.immunity;
//        }
//        else if(this.sortBy == Chromosome.BY_RO){ //good to bad - more to less
//            if(c.tempRo - this.tempRo>0){
//                return 1;
//            }else if(c.tempRo - this.tempRo<0){
//                return -1;
//            }else{
//                return 0;
//            }
//        }
//        else //can also be used for "preferenced fitness_ values" - not implemented yet.
//            throw new UnsupportedOperationException("Not supported yet.");
//        
//    }
//
//    @Override
//    public Object clone() {
//        try {
//            Chromosome chromosome = (Chromosome) super.clone();
//            chromosome.externalData_ = (ExternalData) this.externalData_.clone();
//            chromosome.userInput_ = (UserInput) this.userInput_.clone();
//
//            chromosome.vals_ = (ArrayList<Double>)vals_.clone();
//            Collections.copy(chromosome.vals_, vals_);
//
//            chromosome.fitness_ = (ArrayList<Double>)fitness_.clone();
//            Collections.copy(chromosome.fitness_, fitness_);
//
//
//            chromosome.violations_ = (ArrayList<Integer>)violations_.clone();
//            Collections.copy(chromosome.violations_, violations_);
//
//            chromosome.satisfactions_ = (ArrayList<Integer>)satisfactions_.clone();
//            Collections.copy(chromosome.satisfactions_, satisfactions_);
//
//            return chromosome;
//
//        } catch (CloneNotSupportedException e) {
//            return null;
//        }
//    }
//
//
//
//}
//
//
//
//
//
//
//
//
//
//
//
//
///// for resetObjFunction
////<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     
//                    //Time tabling problme
//                    //<<
//            //            ArrayList<Double> satisfyCheck;
//            //            int count;
//            //            satisfyCheck = vals_;
//            //
//            //            for (Double i : vals_) {
//            //                count = 0;
//            //                for (Double j : satisfyCheck) {
//            //                    if(externalData_.isViolated(i.intValue(), j.intValue())){
//            //                        satisfyCheck = new ArrayList<Double>(satisfyCheck.subList(0, count));
//            //                        break;
//            //                    }else{
//            //                        ;
//            //                    }
//            //                    count++;
//            //                }
//            //            }
//            //            vals_ = satisfyCheck;
//                    //>>
//
//                    //MAP COLORING PROBLEM
//                    //<<
//            //            updateObjectiveFunctionVars();
//                    //>>
//
//                    //NQueen Problem
//                    //<<
//            //        this.fitness_.clear();
//            //        this.fitness_.add(Double.NaN);
//            //        this.sortBy = Chromosome.BY_SATISFACTIONS;
//            //        this.satisfactions_.clear();
//            //
//            //        repairVal();
//            //
//            //        for (int i = 0; i < vals_.size(); i++) {
//            //            for (int j = i+1; j < vals_.size(); j++) {
//            //                 if(i!=j){
//            //                    if(externalData_.isViolated(vals_.get(i).intValue(), vals_.get(j).intValue(),Math.abs(i-j))){
//            //                        vals_.remove(j);
//            //                        j--;
//            //                        
//            //                    }else{
//            //                        ;
//            //                    }
//            //                }
//            //            }
//            //        }
//            //
//            //        for (Double v : vals_) {
//            //            satisfactions_.add(v.intValue());
//            //        }
//                    //>>
//
//
//                    //<<My CIRCLEs Test Problem
//            //        this.fitness_ = new ArrayList<Double>();
//            //        this.sortBy = Chromosome.BY_VIOLATIONS; //its a MUST
//            //
//            //        if(this.vals_.size() != this.userInput_.totalDecisionVars){
//            //            return;
//            //        }
//            //
//            //        this.fitness_.add( Math.sqrt(4- Math.pow((Double)this.vals_.get(0)-2,2))+2 );
//            //        this.fitness_.add( -Math.sqrt(4- Math.pow((Double)this.vals_.get(0)-2,2))+2 );
//            //
//            //        this.fitness_.add( Math.sqrt(1- Math.pow((Double)this.vals_.get(0)-4.98,2)) +2);
//            //        this.fitness_.add( -Math.sqrt(1- Math.pow((Double)this.vals_.get(0)-4.98,2)) +2);
//            //
//            //        this.fitness_.add( Math.sqrt(1- Math.pow((Double)this.vals_.get(0)-4.5,2)) +1);
//            //        this.fitness_.add( -Math.sqrt(1- Math.pow((Double)this.vals_.get(0)-4.5,2)) +1);
//            //
//            //        ArrayList<Double> x = this.vals_;
//            //        ArrayList<Double> f = this.fitness_;
//            //
//            //        this.violations_.clear();
//            //
//            //        if((x.get(1) > f.get(0) || x.get(1) < f.get(1) || x.get(0)>4)      ||    (x.get(1) < f.get(1) || x.get(1) > f.get(0) || x.get(0)>4 || x.get(0)<0))
//            //            this.violations_.add(1);
//            //
//            //        if((x.get(1) > f.get(2) || x.get(1) < f.get(3) || x.get(0)<3.98)    ||    (x.get(1) < f.get(3) || x.get(1) > f.get(2) || x.get(0)<3.98 || x.get(0)>5.98))
//            //            this.violations_.add(2);
//            //
//            //        if((x.get(1) > f.get(4) || x.get(1) < f.get(5) || x.get(0)<3.5)    ||    (x.get(1) < f.get(5) || x.get(1) > f.get(4) || x.get(0)<3.5 || x.get(0)>5.5))
//            //            this.violations_.add(3);
//                    //>>
//            
////resetObjFunction
////>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
//
//
//
////for updateObjFunction<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//
//                //resetObjectiveFunctionVars();
//
//                //<< Nqueen problem
//        //        this.fitness_.clear();
//        //        this.fitness_.add(Double.NaN);
//        //        this.sortBy = Chromosome.BY_SATISFACTIONS;
//        //        this.satisfactions_.clear();
//        //
//        //        repairVal();
//        //
//        //        int i;
//        //        for (int j = vals_.size()-1; j >=0; j--) {
//        //            i = vals_.size()-1;
//        //            if(i!=j){
//        //                if(externalData_.isViolated(vals_.get(i).intValue(), vals_.get(j).intValue(), i-j)){
//        //                    vals_.remove(i);
//        //                    break;
//        //                }
//        //            }
//        //        }
//        //
//        //
//        //        for (Double v : vals_) {
//        //            satisfactions_.add(v.intValue());
//        //        }                    
//                //>>
//
//
//            //<< NqueenII            
//        //        this.fitness_.clear();
//        //        this.fitness_.add(Double.NaN);
//        //        this.sortBy = Chromosome.BY_SATISFACTIONS;
//        //        this.satisfactions_.clear();
//        //
//        //        repairVal();
//        //
//        //        for (int i = 0; i < vals_.size(); i++) {
//        //            for (int j = i+1; j < vals_.size(); j++) {
//        //                 if(i!=j){
//        //                    if(externalData_.isViolated(vals_.get(i).intValue(), vals_.get(j).intValue(),Math.abs(i-j))){
//        //                        vals_.remove(j);
//        //                        j--;
//        //                    }else{
//        //                        ;
//        //                    }
//        //                }
//        //            }
//        //        }
//        //
//        //        for (Double v : vals_) {
//        //            satisfactions_.add(v.intValue());
//        //        }                    
//            //>>
//
//
//                //<< Graph coloring problem
//        //        this.fitness_.clear();
//        //        this.fitness_.add(Double.NaN);
//        //        this.sortBy = Chromosome.BY_SATISFACTIONS;
//        //        this.satisfactions_.clear();
//        //
//        //        repairVal();
//        //
//        //        for (int i = 0; i < vals_.size(); i++) {
//        //            for (int j = 0; j < vals_.size(); j++) {
//        //                 if(i!=j){
//        //                    if(externalData_.isViolated(vals_.get(i), vals_.get(j))){
//        //                        vals_.remove(j);
//        //                        j--;
//        //                    }else{
//        //                        ;
//        //                    }
//        //                }
//        //            }
//        //        }
//        //
//        //        for (Double v : vals_) {
//        //            satisfactions_.add(v.intValue());
//        //        }
//                //>>
//
//
//                //<<Time Tabling Problem
//        //        this.fitness_.clear();
//        //        this.fitness_.add(Double.NaN); // fitness is not needed.
//        //        this.sortBy = Chromosome.BY_SATISFACTIONS;
//        //
//        //        Integer element;
//        //        boolean bviolated;
//        //
//        //        try {
//        //            if(satisfactions_.isEmpty()){
//        //                if(vals_.size() == 1){
//        //                    satisfactions_.add(vals_.get(0).intValue());
//        //                    return;
//        //                }else{
//        //                    throw new Exception("Incompatible values. Check satisfaction value and vals...");
//        //                }
//        //            }else{ //Only check the last element ..... it will be fast
//        //                element = vals_.get(vals_.size()-1).intValue();//last element
//        //
//        //                bviolated = false;
//        //                for (int i = 0; i < satisfactions_.size(); i++) {
//        //                    if(externalData_.isViolated(satisfactions_.get(i).intValue(), element)){
//        //                        bviolated = true;
//        //                        break;
//        //                    }
//        //                }
//        //
//        //                if(!bviolated){
//        //                    //shift the element from last position to first position
//        //                    this.vals_.add(0,element.doubleValue());//insert to first -NOTE - order doesn't matter in this particular problem
//        //                    this.vals_.remove(this.vals_.size()-1);//then remove from last
//        //                    this.satisfactions_.add(0, element);
//        //                }
//        //            }
//        //        } catch (Exception e) {
//        //            e.printStackTrace();
//        //            Application.getInstance().exit();
//        //        }
//
//
//                //>>
//
//
//
////updateObjFunction>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
//
//
