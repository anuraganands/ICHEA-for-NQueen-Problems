/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package csp;

import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author s425770
 */


class Constraint{
    public double val;
    public boolean isViolated;
    public Constraint(double val){
        this.val = val;
        isViolated = false; //default                 
    }
}

public class CCSPfns {
    static ArrayList<Double> vals;
    static ArrayList<Double> fitness_;
    static ArrayList<Integer> violation;
    static UserInput userInput_;
    static double maxCSPval;
    static ArrayList<Integer> position;
    static ArrayList<ArrayList<ArrayList<Double>>> CSPsols;
    static ArrayList<ArrayList<Double>> tabuConstraints = new ArrayList<ArrayList<Double>>();; //additional dynamic constraints.
    static int maxCSPsolsSize;
    static int lastIdx;
    static ArrayList<Constraint> consFn;
    static double curFitness;
    static int MAX_FUNCTIONAL_CONSTRAINTS;
    
    private CCSPfns(){
        ;
    }
    
    public static void objFn(int fn, final ArrayList<Double> lvals, final ArrayList<Double> lfitness_, 
            final ArrayList<Integer> lconstraints, UserInput luseInput, double lmaxCSPval)
            throws SolutionFoundException{ 
//            ArrayList<ArrayList<ArrayList<Double>>> lCSPsols){
        vals = lvals;
        fitness_ = lfitness_;
        violation = lconstraints;
        userInput_ = luseInput;
        maxCSPval = lmaxCSPval;
        position = new ArrayList<Integer>(); //lposition;  
        CSPsols= CspProcess.CSPsols;
        maxCSPsolsSize = userInput_.population; //(int)Math.ceil(1.0*userInput_.population/userInput_.totalConstraints);
        consFn = new ArrayList<Constraint>();
        lastIdx = getLastIdx(); //min should be 1 as defined in CspProcess initialize() method.
        curFitness = Double.NaN;
        MAX_FUNCTIONAL_CONSTRAINTS = userInput_.totalConstraints - userInput_.totalDecisionVars; // becuase I take total_constraint as total decision vars + given constraints
//        tabuConstraints = new ArrayList<ArrayList<Double>>();
        
        if(isCSPsolsAvailable(CSPsols.get(lastIdx)) && !CSPsols.get(lastIdx).isEmpty()){
            CSPsols.add(new ArrayList<ArrayList<Double>>());
            for (int i = 0; i < userInput_.totalConstraints; i++) {
                CSPsols.get(lastIdx+1).add(new ArrayList<Double>());
            }
        }
        lastIdx = getLastIdx();  
        
        violation.clear();
        
        switch (fn){
            case 1:
                broyden10(CspProcess.dynamicConstraintNo);//hs109;// h77 //G5(); //neuro(); //eco9(); //G5(); //broyden10();//chem(); // //TP7();// 
                break;
            case 2:
                //?
                break;
            case 3:
                //?
                break;
            case 4:
                //?
                break;
            default:
                throw new UnsupportedOperationException("Function not found!");
        }
        
                fitness_.clear();
        
        double fit=0;
        if(userInput_.bWeighted){
            if(maxCSPval>= Double.MAX_VALUE)
                fitness_.add(violation.size()*1.0);
            else
                fitness_.add(curFitness+violation.size()*Math.abs(maxCSPval));
        }else{
            fit = 0;
            for (Constraint cns : consFn) {
                if(cns.isViolated){//violation
                    fit+= Math.abs(cns.val);
                }
            }
            fitness_.add(fit);
            
//            fitness_.add(violation.size()*1.0);
        }
        
        fitness_.add(curFitness); //fitness_[1]   
        
        
        //csp sol
        if(violation.size() == 0 && CspProcess.dynamicConstraintNo >= MAX_FUNCTIONAL_CONSTRAINTS && !CspProcess.bInTransition){            
            throw new SolutionFoundException("At least one CSP solution found");
        }
    }
    
    /**
     * This one adds dynamic constraints to avoid system being stuck into local optima
     * Addition a of new dynamic constraint will keep the system away from revisiting it.
     * 
     * @param center center of your defined hyper-sphere
     * @param radus radius of hyper-sphere
     */
    public static void addTabuConstraint(ArrayList<Double> center, double radius){
                
        // get the middle point between bestSoFarVal and nearestOfBestSoFar
        //center = MyMath.vectorAddition(bestSoFarVal, nearestOfBestSoFar);
        //center = MyMath.constMultiplicationToVector(0.5, center); 
        
        
        //creat a hyper-sphere where contraint is inner points, hence outer points are allowed
        //(x1-c1)^2 + (x2-c2)^2 + .... +(xn-cn)^2 - r^2 >= 0
        //Store constants c1,c2...,ci,...,cn which is a centre, ofcourse
        tabuConstraints.add(new ArrayList<Double>());
        
        for (Double ci : center) {
            tabuConstraints.get(tabuConstraints.size()-1).add(ci);
        }
        tabuConstraints.get(tabuConstraints.size()-1).add(radius); //last one is radius
        
        
    }
    
    private static boolean isCSPsolsAvailable(ArrayList<ArrayList<Double>> oneCSPsols){
        boolean avail = true;
        
        for (ArrayList<Double> valGrp: oneCSPsols) {
            if(valGrp.isEmpty()){
                avail = false;
                break;
            }
        }
        
        return avail;
    }
    
    /**
     * If the last idx is full the first one is dequeued to make space for the new one.
     * @return 
     */
    private static int getLastIdx(){
        int lastIdx = CSPsols.size()-1;
        lastIdx = lastIdx % maxCSPsolsSize;
        
        if(lastIdx == maxCSPsolsSize-1){ //the last one
            CSPsols.remove(0);
            lastIdx = CSPsols.size()-1;
            lastIdx = lastIdx % maxCSPsolsSize;
        }
        
        return lastIdx;
    }
    
    /**
     * function from : Evolutionary Algorithms for Constrained Parameter Optimization Problems
     * Author: Zbigniew Michalewicz and Marc Schoenauery
     * function name G1
     */
    private static void G1(){
        Double x[] = new Double[userInput_.totalDecisionVars];        
        vals.toArray(x);
        // total vars = 13
        // total constraints = 9
        //0 <= xi <= 1, i = 0-8, 
        //0 <= xi <= 100, i = 9,10,11
        //0 <= x12 <= 1.
        
        
        
        double sum1, sum2;
        
        sum1 = 0;
        for (int i = 0; i <= 3; i++) {
            sum1 += x[i]*x[i];
        }
        sum2 = 0;
        for (int i = 4; i <= 12; i++) {
            sum2 += x[i];
        }
        
        curFitness = 5*x[0] + 5*x[1] + 5*x[2] + 5*x[3] -5*sum1 - sum2;
        violation.clear();       
        
        final int totalCons = 9;          
        double t[] = new double[totalCons];

        t[0] = -(2*x[0] + 2*x[1] + x[9] + x[10] -10); //>=0
        t[1] = -(2*x[0] + 2*x[2] + x[9] + x[11] -1);//>=0
        t[2] = -(2*x[1] + 2*x[2] + x[10] + x[11] -10);//>=0
        t[3] = -(-8*x[0] + x[9]);//>=0
        t[4] = -(-8*x[1] + x[10]);//>=0
        t[5] = -(-8*x[2] + x[11]);//>=0
        t[6] = -(-2*x[3] - x[4] + x[9]);//>=0
        t[7] = -(-2*x[5] - x[6] + x[10]);//>=0
        t[8] = -(-2*x[7] - x[8] + x[11]);//>=0
        
      
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
                
        
        for (int i = 0; i < t.length; i++) {
            if(consFn.get(i).val >= 0.0){
//                if(CSPsols.size() != maxCSPsolsSize-1 && CSPsols.get(lastIdx).get(i).isEmpty())
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
 
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = t.length; i < t.length+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }                      
    }
    
    private static void G2(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 2 (prefered)
        // total constraints = 2+2
        //0 <= xi <= 10, i = 0-n, 

        
        double sum1, product1, sum2;
        final int n = userInput_.totalDecisionVars;
        
        sum1 = 0.0;
        for (int i = 0; i < n; i++) {
            sum1 += Math.pow(Math.cos(x[i]),4);
        }
        
        
        product1 = 1.0;
        for (int i = 0; i < n; i++) {
            product1 *= Math.pow(Math.cos(x[i]),2);
        }
        
        sum2 = 0;
        for (int i = 0; i < n; i++) {
            sum2 += (i+1)*x[i]*x[i];
        }
        
        curFitness = Math.abs((sum1-2*product1)/Math.sqrt(sum2));
        
        
        violation.clear();
        
        
//        x[0] = 7.872373; x[1] = 7.408733; x[2] = 3.150534; x[3] = 0.384384; x[4] = 2.187453;
//x[5] = 7.330112; x[6] = 9.641181; x[7] = 0.631615; x[8] = 5.632460; x[9] = 8.265453;
//x[10] = 9.135217; x[11] = 6.150099; x[12] = 3.438410; x[13] = 4.274436; x[14] = 6.636882;
//x[15] = 3.490745; x[16] = 5.703963; x[17] = 9.155801; x[18] = 2.181950; x[19] = 2.930025;
        
        final int totalCons = 2;          
        double t[] = new double[totalCons];

        product1 = 1.0;
        for (int i = 0; i < n; i++) {
            product1 *= x[i];
        }
        
        sum1 = 0.0;
        for (int i = 0; i < n; i++) {
            sum1 += x[i];
        }        
        
        t[0] = product1-0.75; //>=0
        t[1] = 7.5*n - sum1;//>=0
        
      
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
                
        
        for (int i = 0; i < t.length; i++) {
            if(consFn.get(i).val >= 0.0){
//                if(CSPsols.size() != maxCSPsolsSize-1 && CSPsols.get(lastIdx).get(i).isEmpty())
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
 
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = t.length; i < t.length+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }                      
              
    }
    
    private static void G10(){
         Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 8
        // total constraints = 6+8
        // 100<=xi<=10,000 i = 1
        // 1000<=xi<=10,000 i = 2,3
        // 10<=xi<=1000 i = 4..8
        
        curFitness = 0.0;
        violation.clear();
        double t[] = new double[6];
        t[0] = -1.0 + (0.0025*(x[3]+x[5])); //<=0
        t[1] = -1.0 + (0.0025*(x[4]+x[6]-x[3]));//<=0
        t[2] = -1.0 + (0.01*(x[7]-x[4]));//<=0
        t[3] = (-x[0]*x[5])+(833.33252*x[3])+(100.0*x[0])-(83333.333);//<=0
        t[4] = (-x[1]*x[6])+(1250.0*x[4])+(x[1]*x[3])-(1250.0*x[3]);//<=0
        t[5] = (-x[2]*x[7])+(1250000)+(x[2]*x[4])-(2500.0*x[4]);//<=0
        
        t[0] = -t[0];
        t[1] = -t[1];
        t[2] = -t[2];
        t[3] = -t[3];
        t[4] = -t[4];
        t[5] = -t[5];
        
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
        
        for (int i = 0; i < consFn.size(); i++) {
            if(consFn.get(i).val >= 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
         
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
      
        for (int i = t.length; i < t.length+x.length; i++) {
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        } 
    }
    
    private static void G9(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 7
        // total constraints = 4 + 7
        // -10<=xi<=10, i = 0..6
        
        curFitness = (x[0]-10)*(x[0]-10) + 5*(x[1]-12)*(x[1]-12) + Math.pow(x[2], 4)
                +3*(x[3]-11)*(x[3]-11) +10*Math.pow(x[4], 6)+ 7*x[5]*x[5] + Math.pow(x[6], 4)
                -4*x[5]*x[6] - 10*x[5] - 8*x[6];
        violation.clear();
        
        double t[] = new double[4];
        t[0] = 127-2*x[0]*x[0] - 3*Math.pow(x[1], 4) -x[2] - 4*Math.pow(x[3], 2)-5*x[4];//>=0;
        t[1] = 282-7*x[0]-3*x[1]-10*x[2]*x[2]-x[3]+x[4];//>=0;
        t[2] = 196-23*x[0]-x[1]*x[1]-6*x[5]*x[5]+8*x[6];//>=0;
        t[3] = -4*x[0]*x[0] - x[1]*x[1]+3*x[0]*x[1] - 2*x[2]*x[2] - 5*x[5] + 11*x[6];// >=0;
        
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
        
        for (int i = 0; i < consFn.size(); i++) {
            if(consFn.get(i).val >= 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
         
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
      
        for (int i = t.length; i < t.length+x.length; i++) {
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        } 
        

    }
    
    
    private static void G7(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 10
        // total constraints = 8+10
        // -10.0 <=xi <= 10.0, i = 0..9
        
        curFitness = x[0]*x[0]+ x[1]*x[1]+ x[0]*x[1] - 14*x[0] - 16*x[1] + (x[2]-10)*(x[2]-10) + 
                4*(x[3] - 5)*(x[3] - 5) + (x[4]- 3)*(x[4]-3) + 2*(x[5]-1)*(x[5]-1) +5*x[6]*x[6]+ 
                7*(x[7]- 11)*(x[7]- 11) + 2*(x[8]- 10)*(x[8]- 10) + (x[9]- 7)*(x[9]- 7) + 45;
        
        double t[] = new double[8];
        t[0] = 105 - 4*x[0] - 5*x[1] + 3*x[6] - 9*x[7];// >= 0
        t[1] = -3*(x[0]-2)*(x[0]-2) - 4*(x[1] - 3)*(x[1] - 3) - 2*x[2]*x[2]+ 7*x[3] + 120;// >=0
        t[2] = -10*x[0] + 8*x[1] + 17*x[6] - 2*x[7];// >=0
        t[3] = -x[0]*x[0] - 2*(x[1]- 2)*(x[1]- 2) + 2*x[0]*x[1] - 14*x[4] + 6*x[5];// >= 0,
        t[4] = 8*x[0] - 2*x[1] - 5*x[8] + 2*x[9] + 12;// >= 0
        t[5] = -5*x[0]*x[0] - 8*x[1] - (x[2] - 6)*(x[2] - 6) + 2*x[3] + 40;// >= 0,
        t[6] = 3*x[0] - 6*x[1] - 12*(x[8] - 8)*(x[8] - 8) + 7*x[9];// >= 0
        t[7] = -0.5*(x[0] - 8)*(x[0] - 8) - 2*(x[1] - 4)*(x[1] - 4) - 3*x[4]*x[4]+ x[5] + 30;// >= 0
        
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
        
        
        //all >= 0
        for (int i = 0; i < consFn.size(); i++) {
             if(consFn.get(i).val >= 0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            }
        }
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
      
        for (int i = t.length; i < t.length+x.length; i++) {
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }
    }
    
    
    private static void G5(final int tm) throws SolutionFoundException{
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 4
        // total constraints = 5+4
        // 0 <= xi <= 1200, i = 0, 1
        // -0.55<=xi<=0.55, i = 2, 3
        
        curFitness = 3*x[0] + 0.000001*x[0]*x[0]*x[0] + 2*x[1] + 0.000002/3*x[1]*x[1]*x[1];
        violation.clear();
        
        
        final int minPow = -5;
        final int maxPow = 2;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
//        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
//        for (int i = 1; i < total_divisions; i++) {
//            if(lastIdx>grpLenght*i){
//                grp++;
//                continue;
//            }else{
//                break;
//            }
//        }
        grp = CspProcess.negFeasibleRange;//Math.min(grp,CspProcess.negFeasibleRange);
               
        final double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            CspProcess.bInTransition = false;
        }
        
        final int tabuSize = tabuConstraints.size();
        final int totalFunctionalConstraints = 5;        
        double t[] = new double[totalFunctionalConstraints+tabuSize];
        
        t[0] = x[3]-x[2]+0.55;//>=0;
        t[1] = x[2]-x[3]+0.55;//>=0;
        t[2] = 1000*Math.sin(-x[2]-0.25)+1000*Math.sin(-x[3]-0.25)+894.8-x[0];//=0;
        t[3] = 1000*Math.sin(x[2]-0.25)+1000*Math.sin(x[2]-x[3]-0.25)+894.8-x[1];//=0;
        t[4] = 1000*Math.sin(x[3]-0.25)+1000*Math.sin(x[3]-x[2]-0.25)+1294.8;//=0;        
        
        //Short term memory - for more than 5 tabus start deleting/forgetting the oldest one.
        final int SHORT_TERM_MEMORY = 5;
        for (int i = SHORT_TERM_MEMORY; i < tabuSize; i++) {
            t[totalFunctionalConstraints+(i-SHORT_TERM_MEMORY)] = 1.0; //to make it >=0
        }
        
        double sum=0.0;
        if(tabuSize>0){
            for (int cons = 0; cons < tabuConstraints.size(); cons++) {
                ArrayList<Double> centers = tabuConstraints.get(cons);
                if(x.length != centers.size()-1){ //last one is radius
                    throw new SolutionFoundException("NO SOLUTION! BUG!!!! Center and varX size must be same. Check your code");
                }
                sum = 0.0;
                //This is equation of a hyper-sphere
                for (int i = 0; i < x.length; i++) {
                    sum += Math.pow(x[i]-centers.get(i),2);                                        
                }
                //hyper-sphere equation. (x1-c1)^2 + (x2-c2)^2 + .... +(xn-cn)^2 - r^2 >= 0
                t[totalFunctionalConstraints+cons] = sum - Math.pow(centers.get(centers.size()-1),2); //>=0
            }
        }
            
        final int USED_CONS = Math.min(tm,t.length);

      
        int consOrder[] = new int[100];//{0,1,2,3,4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, };//{0,3,4,5,6,9,10,8,1,2,7};
        Arrays.fill(consOrder, -1);
        int consOrderFncons[] = {0,1,2,3,4};
        for (int i = 0; i < consOrderFncons.length; i++) {
            consOrder[i] = consOrderFncons[i];
            
        }        
        for (int i = totalFunctionalConstraints; i < t.length; i++) {
            consOrder[i] = i;            
        }
//        int consOrder[] = {0};
             
        for (int i = 0; i < USED_CONS; i++) {
            consFn.add(new Constraint(t[consOrder[i]]));            
        }
        

        int idx; 
        for (int i = 0; i < USED_CONS; i++) {
            idx = consOrder[i];
          
            if((idx>=0 && idx<= 1) || idx >= totalFunctionalConstraints){
                if(consFn.get(i).val >= 0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                }
            }else{
                if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                } 
            }
        }
        
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = USED_CONS; i < USED_CONS+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-USED_CONS) || 
            consFn.get(i).val>userInput_.maxVals.get(i-USED_CONS)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }          
    }
    
    private static void G4(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 5
        // total constraints = 3+5
        //78 <= x0 <= 102
        //33 <= x1 <= 45
        //27 <= xi <= 45 ,i = 2, 3, 4      
       
    
        curFitness = 5.3578547*x[2]*x[2] + 0.8356891*x[0]*x[4] + 37.293239*x[0] - 40792.141; //need to maximize        
        
        violation.clear();
        
        double t[] = new double[3];
        t[0] = 85.334407 + 0.0056858*x[1]*x[4] + 0.00026 *x[0]*x[3] - 0.0022053*x[2]*x[4];
        t[1] = 80.51249 + 0.0071317*x[1]*x[4] + 0.0029955*x[0]*x[1] + 0.0021813*x[2]*x[2];
        t[2] = 9.300961 + 0.0047026*x[2]*x[4] + 0.0012547*x[0]*x[2] + 0.0019085*x[2]*x[3];        
        
        consFn.add(new Constraint(t[0])); //==0
        consFn.add(new Constraint(t[1])); //==0
        consFn.add(new Constraint(t[2])); //==0
        
        if( 0.0 <= consFn.get(0).val && consFn.get(0).val <=92.0){
            CSPsols.get(lastIdx).set(0, (ArrayList<Double>)vals.clone());            
        }else{
            consFn.get(0).isViolated = true;
            violation.add(0);
        }                  
        
        if( 90.0 <= consFn.get(1).val && consFn.get(1).val <=110.0){
            CSPsols.get(lastIdx).set(1, (ArrayList<Double>)vals.clone());
        }else{
            consFn.get(1).isViolated = true;
            violation.add(1);            
        } 
                
        if( 20.0 <= consFn.get(2).val && consFn.get(2).val <=25.0){
            CSPsols.get(lastIdx).set(2, (ArrayList<Double>)vals.clone());
        }else{
            consFn.get(2).isViolated = true;
            violation.add(2);            
        } 
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = t.length; i < t.length+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }     
        
        
    }
    
    private static void G8(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 2 
        // total constraints = 2+2
        // 0 <= x0 <= 10,
        // 0 <= x1 <= 10,
        
        curFitness = 0.0; 
        violation.clear();
        
        double t[] = new double[2];
        t[0] = 1.0-x[0]*x[0] + x[1]*x[1];//>=0
        t[1] = x[0] -(x[1]-4)*(x[1]-4)-1;//>=0
       
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
        
        for (int i = 0; i < consFn.size(); i++) {
            if(consFn.get(i).val >= 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
         
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
      
        for (int i = t.length; i < t.length+x.length; i++) {
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }
    }
    
    
    private static void G6(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 2 
        // total constraints = 2+2
        // 13 <= x0 <= 100,
        // 0 <= x1 <= 100,
        
        curFitness = Math.pow(x[0]-10,3) + Math.pow(x[1] - 20,3);
        violation.clear();
        
        double t[] = new double[2];
        t[0] = Math.pow(x[0]-5,2)+Math.pow(x[1]-5,2)-100;//>=0;
        t[1] = -Math.pow(x[0]-6,2)-Math.pow(x[1]-5,2)+82.81;//>=0;
       
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
        
        for (int i = 0; i < consFn.size(); i++) {
            if(consFn.get(i).val >= 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
      
        for (int i = t.length; i < t.length+x.length; i++) {
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }
        
    }
    
     private static void G3(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 10 (prefered)
        // total constraints = 1+10
        //0 <= xi <= 1, i = 0-n, 

        
        double sum1, product1;
        int n = userInput_.totalDecisionVars;
        
        product1 = 1;
        for (int i = 0; i < n; i++) {
            product1 *= x[i];
        }
    
        curFitness = -Math.pow(Math.sqrt(n),n)*product1; //need to maximize                
        violation.clear();
        
        final int minPow = -3;
        final int maxPow = 2;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
        for (int i = 1; i < total_divisions; i++) {
            if(lastIdx>grpLenght*i){
                grp++;
                continue;
            }else{
                break;
            }
        }
             
        
        final double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            CspProcess.bInTransition = false;
        }

        
        final int totalCons = 1;
        double t[] = new double[totalCons];
        
        
        sum1 = 0;
        for (int i = 0; i < n; i++) {
            sum1 += x[i]*x[i];
        }
        
        t[0] = sum1-1;
        
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
                
        for (int i = 0; i < consFn.size(); i++) {
            if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = t.length; i < t.length+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }                     
    }
    
    
    
     private static void distributedCircles(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 2 
        // total constraints = 5
        //-110 <= xi <= 110, i = 0-n, 

        
        double sum1, product1;
        int n = userInput_.totalDecisionVars;
        double r = 2.0;
        
        consFn.add(new Constraint((x[0]-97)*(x[0]-97) + x[1]*x[1]-r*r));// cons 1 <=0
	consFn.add(new Constraint((x[0]-100)*(x[0]-100) + x[1]*x[1]-r*r));//cons 2 <=0
 
	consFn.add(new Constraint((x[0]+95)*(x[0]+95) + x[1]*x[1]-r*r));//cons 1 <=0
	consFn.add(new Constraint((x[0]+100)*(x[0]+100) + x[1]*x[1]-r*r));//cons 1 <=0
 
	consFn.add(new Constraint((x[0]+50)*(x[0]+50) + (x[1]-50)*(x[1]-50)-r*r));//cons 2 <=0
                          
    
        curFitness = 0.0; //need to maximize
        
        
        violation.clear();
        
        if( consFn.get(0).val<=0 || consFn.get(2).val<=0 || consFn.get(3).val<=0){
            CSPsols.get(lastIdx).set(0, (ArrayList<Double>)vals.clone());
        }else{
            consFn.get(0).isViolated = true;
            violation.add(0);            
        }   
        

        if( consFn.get(1).val>0 || consFn.get(4).val>0){
            CSPsols.get(lastIdx).set(0, (ArrayList<Double>)vals.clone());
        }else{
            consFn.get(0).isViolated = true;
            violation.add(1);            
        } 
                               
    }
     
    private static void TP7(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 5
        // total constraints = 3
        // -2.3 <=xi <= 2.3, i = 0, 1   
        // -3.2 <=xi <= 3.2, i = 2,3,4 
        //-2.3, -2.3, -3.2, -3.2, -3.2
        //2.3, 2.3, 3.2, 3.2, 3.2
        
        curFitness = Math.exp(x[0]*x[1]*x[2]*x[3]*x[4]);
        violation.clear();
        
        double t[] = new double[3];
        t[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4]-10;//=0
        t[1] = x[1]*x[2] -5*x[3]*x[4];//=0;
        t[2] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1] +1 ;//=0;

        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
        
        for (int i = 0; i < consFn.size(); i++) {
            if(Math.abs(MyMath.roundN(consFn.get(i).val,2)) == 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
        
    } 
    
    private static void broyden10(final int tm) throws SolutionFoundException{
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 10
        // total constraints = 10+10
		
        // -1.0E8<= x[i] <= 1.0E8; i = 0..9 
 
        curFitness = 0.0;
        violation.clear();
        
        final int minPow = -1;
        final int maxPow = 20;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
//        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
//        for (int i = 1; i < total_divisions; i++) {
//            if(lastIdx>grpLenght*i){
//                grp++;
//                continue;
//            }else{
//                break;
//            }
//        }
            
        grp = CspProcess.negFeasibleRange;//Math.min(grp,CspProcess.negFeasibleRange);
        
        
        double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            dt = MIN_DT;
            CspProcess.bInTransition = false;
        }
        
        final int tabuSize = tabuConstraints.size();
        final int totalFunctionalConstraints = 10;
        
        double t[] = new double[totalFunctionalConstraints+tabuSize];
//        double t[] = new double[1];
        
        t[0] = x[0]*(2+5*x[0]*x[0]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1])); //= 0
        t[1] = x[1]*(2+5*x[1]*x[1]) + 1 - (x[0]*(1+x[0]) + x[2]*(1+x[2])); //= 0
        t[2] = x[2]*(2+5*x[2]*x[2]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[3]*(1+x[3])); //= 0
        t[3] = x[3]*(2+5*x[3]*x[3]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[4]*(1+x[4])); //= 0
        t[4] = x[4]*(2+5*x[4]*x[4]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[3]*(1+x[3])+ x[5]*(1+x[5])); //= 0

        t[5] = x[5]*(2+5*x[5]*x[5]) + 1 - (x[0]*(1+x[0]) + x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[3]*(1+x[3])+ x[4]*(1+x[4]) + x[6]*(1+x[6])) ; //= 0
        t[6] = x[6]*(2+5*x[6]*x[6]) + 1 - (x[1]*(1+x[1]) + x[2]*(1+x[2]) + x[3]*(1+x[3]) + x[4]*(1+x[4])+ x[5]*(1+x[5]) + x[7]*(1+x[7])) ; //= 0
        t[7] = x[7]*(2+5*x[7]*x[7]) + 1 - (x[2]*(1+x[2]) + x[3]*(1+x[3]) + x[4]*(1+x[4]) + x[5]*(1+x[5])+ x[6]*(1+x[6]) + x[8]*(1+x[8])) ; //= 0
        t[8] = x[8]*(2+5*x[8]*x[8]) + 1 - (x[3]*(1+x[3]) + x[4]*(1+x[4]) + x[5]*(1+x[5]) + x[6]*(1+x[6])+ x[7]*(1+x[7]) + x[9]*(1+x[9])) ; //= 0
        t[9] = x[9]*(2+5*x[9]*x[9]) + 1 - (x[4]*(1+x[4]) + x[5]*(1+x[5]) + x[6]*(1+x[6]) + x[7]*(1+x[7])+ x[8]*(1+x[8])) ; //= 0
    
        //Short term memory - for more than 5 tabus start deleting/forgetting the oldest one.
        final int SHORT_TERM_MEMORY = 5;
        for (int i = SHORT_TERM_MEMORY; i < tabuSize; i++) {
            t[totalFunctionalConstraints+(i-SHORT_TERM_MEMORY)] = 1.0; //to make it >=0
        }
        
        double sum=0.0;
        if(tabuSize>0){
            for (int cons = 0; cons < tabuConstraints.size(); cons++) {
                ArrayList<Double> centers = tabuConstraints.get(cons);
                if(x.length != centers.size()-1){ //last one is radius
                    throw new SolutionFoundException("NO SOLUTION! BUG!!!! Center and varX size must be same. Check your code");
                }
                sum = 0.0;
                //This is equation of a hyper-sphere
                for (int i = 0; i < x.length; i++) {
                    sum += Math.pow(x[i]-centers.get(i),2);                                        
                }
                //hyper-sphere equation. (x1-c1)^2 + (x2-c2)^2 + .... +(xn-cn)^2 - r^2 >= 0
                t[totalFunctionalConstraints+cons] = sum - Math.pow(centers.get(centers.size()-1),2); //>=0
            }
        }
               
        final int USED_CONS = Math.min(tm,t.length);

      
        int consOrder[] = new int[100];//{0,1,2,3,4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, };//{0,3,4,5,6,9,10,8,1,2,7};
        Arrays.fill(consOrder, -1);
        int consOrderFncons[] = {0,1,2,3,4,5,6,7,8,9};
        for (int i = 0; i < consOrderFncons.length; i++) {
            consOrder[i] = consOrderFncons[i];
            
        }        
        for (int i = totalFunctionalConstraints; i < t.length; i++) {
            consOrder[i] = i;            
        }
//        int consOrder[] = {0};
       
        
        for (int i = 0; i < USED_CONS; i++) {
            consFn.add(new Constraint(t[consOrder[i]]));            
        }
        
        
        int idx; 
        for (int i = 0; i < USED_CONS; i++) {
            idx = consOrder[i];
          
            if(idx >= totalFunctionalConstraints){
                if(consFn.get(i).val >= 0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                }
            }else{
                if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                } 
            }
        }
        
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = USED_CONS; i < USED_CONS+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-USED_CONS) || 
            consFn.get(i).val>userInput_.maxVals.get(i-USED_CONS)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }        
    }
    
    
    private static void chem(){
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 5
        // total constraints = 5+5		
        // 0<= x[i] <= 1.0E8; i = 0..4
        
        curFitness = 0.0;
        violation.clear();
        
//        double dt = 1.0E-1;       

        final int minPow = -1;
        final int maxPow = 2;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
        for (int i = 1; i < total_divisions; i++) {
            if(lastIdx>grpLenght*i){
                grp++;
                continue;
            }else{
                break;
            }
        }
             
        if(lastIdx == 0){
            lastIdx = lastIdx;
        }
        
        final double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            CspProcess.bInTransition = false;
        }

        
        
        double r = 10;
        double r5 = 0.193;
        double r6= 0.002597/Math.sqrt(40.0);
        double r7 = 0.003448/Math.sqrt(40.0);
        double r8 = 0.00001799/40;
        double r9 = 0.0002155/Math.sqrt(40.0);
        double r10= 0.00003846/40;

        double t[] = new double[5];

        t[0] = 3*x[4] - x[0]*(x[1] + 1); //=0
        t[1] = x[1]*(2*x[0] + x[2]*x[2] + r8 + 2*r10*x[1] + r7*x[2] + r9*x[3]) + x[0] - r*x[4]; //=0
        t[2] = x[2]*(2*x[1]*x[2] + 2*r5*x[2] + r6 + r7*x[1]) - 8*x[4]; //=0
        t[3] = x[3]*(r9*x[1] + 2*x[3]) - 4*r*x[4]; //=0
        t[4]  = x[1]*(x[0]  + r10*x[1] + x[2]*x[2] + r8  +r7*x[2] + r9*x[3]) + x[0] + r5*x[2]*x[2] + x[3]*x[3] - 1 + r6*x[2];//=0

        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
                
        for(int i = 0; i< t.length; i++){
            if(dt - Math.abs(consFn.get(i).val) >= 0.0){
//                if(CSPsols.size() != maxCSPsolsSize-1 && CSPsols.get(lastIdx).get(i).isEmpty())
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = t.length; i < t.length+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
//                 if(CSPsols.size() != maxCSPsolsSize-1 && CSPsols.get(lastIdx).get(i).isEmpty())
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }
        

    }
    
    
    private static void chem(final int tm) throws SolutionFoundException{
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 5
        // total constraints = 5+5		
        // 0<= x[i] <= 1.0E8; i = 0..4
        
        curFitness = 0.0;
        violation.clear();      

        final int minPow = -1;
        final int maxPow = 2;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
        for (int i = 1; i < total_divisions; i++) {
            if(lastIdx>grpLenght*i){
                grp++;
                continue;
            }else{
                break;
            }
        }
         
        grp = Math.min(grp,CspProcess.negFeasibleRange);
                
        double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            dt = MIN_DT;
            CspProcess.bInTransition = false;
        }

        
        
        double r = 10.0;
        double r5 = 0.193;
        double r6= 0.002597/Math.sqrt(40.0);
        double r7 = 0.003448/Math.sqrt(40.0);
        double r8 = 0.00001799/40;
        double r9 = 0.0002155/Math.sqrt(40.0);
        double r10= 0.00003846/40;

        final int tabuSize = tabuConstraints.size();
        final int totalFunctionalConstraints = 5;        
        double t[] = new double[totalFunctionalConstraints+tabuSize];
//        double t[] = new double[1];

        t[0] = 3*x[4] - x[0]*(x[1] + 1); //=0
        t[1] = x[1]*(2*x[0] + x[2]*x[2] + r8 + 2*r10*x[1] + r7*x[2] + r9*x[3]) + x[0] - r*x[4]; //=0
        t[2] = x[2]*(2*x[1]*x[2] + 2*r5*x[2] + r6 + r7*x[1]) - 8*x[4]; //=0
        t[3] = x[3]*(r9*x[1] + 2*x[3]) - 4*r*x[4]; //=0
        t[4]  = x[1]*(x[0]  + r10*x[1] + x[2]*x[2] + r8  +r7*x[2] + r9*x[3]) + x[0] + r5*x[2]*x[2] + x[3]*x[3] - 1 + r6*x[2];//=0

        double sum=0.0;
        if(tabuSize>0){
            for (int cons = 0; cons < tabuConstraints.size(); cons++) {
                ArrayList<Double> centers = tabuConstraints.get(cons);
                if(x.length != centers.size()-1){ //last one is radius
                    throw new SolutionFoundException("NO SOLUTION! BUG!!!! Center and varX size must be same. Check your code");
                }
                sum = 0.0;
                //This is equation of a hyper-sphere
                for (int i = 0; i < x.length; i++) {
                    sum += Math.pow(x[i]-centers.get(i),2);                                        
                }
                //hyper-sphere equation. (x1-c1)^2 + (x2-c2)^2 + .... +(xn-cn)^2 - r^2 >= 0
                t[totalFunctionalConstraints+cons] = sum - Math.pow(centers.get(centers.size()-1),2); //>=0
            }
        }
        
        //Short term memory - for more than 5 tabus start deleting/forgetting the oldest one.
        final int SHORT_TERM_MEMORY = 5;
        for (int i = SHORT_TERM_MEMORY; i < tabuSize; i++) {
            t[totalFunctionalConstraints+(i-SHORT_TERM_MEMORY)] = 1.0; //to make it >=0
        }
        
        final int USED_CONS = Math.min(tm,t.length);

        int consOrder[] = new int[100];//{0,1,2,3,4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, };//{0,3,4,5,6,9,10,8,1,2,7};
        Arrays.fill(consOrder, -1);
        int consOrderFncons[] = {0,1,2,3,4};
        for (int i = 0; i < consOrderFncons.length; i++) {
            consOrder[i] = consOrderFncons[i];
            
        }        
        for (int i = totalFunctionalConstraints; i < t.length; i++) {
            consOrder[i] = i;            
        }
//        int consOrder[] = {0};
        
        for (int i = 0; i < USED_CONS; i++) {
            consFn.add(new Constraint(t[consOrder[i]]));            
        }
             
        int idx;
        for (int i = 0; i < USED_CONS; i++) {
            idx = consOrder[i];
            
            if(idx>4){ //- see HS109 for this line
                if(consFn.get(i).val >= 0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                }
            }else{
                if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                } 
            }
        }
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = USED_CONS; i < USED_CONS+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-USED_CONS) || 
            consFn.get(i).val>userInput_.maxVals.get(i-USED_CONS)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }
    }
    
    private static void neuro(){
    Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 6
        // total constraints = 9+6
        // -100<= x[i] <= 100; i = all 5

        
        curFitness = 0.0;
        violation.clear();
        
        final int minPow = -3;
        final int maxPow = 2;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
        for (int i = 1; i < total_divisions; i++) {
            if(lastIdx>grpLenght*i){
                grp++;
                continue;
            }else{
                break;
            }
        }
             
        if(lastIdx == 0){
            lastIdx = lastIdx;
        }
        
        final double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            CspProcess.bInTransition = false;
        }

        
        final int totalCons = 9;
        double t[] = new double[totalCons];    
                                
        double C1 = 5.0;
        double C2 = 4.0;
        double C3 = 3.0;
        double C4 = 2.0;
  
        t[0] = x[0]*x[0] + x[2]*x[2]             - 1;//=0
        t[1] = x[1]*x[1] + x[3]*x[3]             - 1;//=0
        t[2] = x[4]*x[0]*x[0]*x[0] + x[5]*x[1]*x[1]*x[1]       - C2;//=0
        t[3] = x[4]*x[0]*x[2]*x[2] + x[5]*x[3]*x[3]*x[1] - C3;//=0
        t[4] = x[4]*x[2]*x[2]*x[2] + x[5]*x[3]*x[3]*x[3]       - C1;//=0
        t[5] = x[4]*x[0]*x[0]*x[2] + x[5]*x[1]*x[1]*x[3] - C4;//=0
        
        t[6] = x[0] - x[1]; //>= 0
        t[7] = x[0]; //>= 0;
        t[8] = x[1]; //>= 0;
        
        
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
                
        for (int i = 0; i < 6; i++) {
            if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
        for (int i = 6; i < t.length; i++) {
            if(consFn.get(i).val >= 0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            }  
        }
        
        
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = t.length; i < t.length+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }
    }
     
    private static void eco9(){
         Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 8
        // total constraints = 8+8
        // -100<= x[i] <= 100; i = all 8

        
        curFitness = 0.0;
        violation.clear();
        
        final int minPow = -1;
        final int maxPow = 2;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
        for (int i = 1; i < total_divisions; i++) {
            if(lastIdx>grpLenght*i){
                grp++;
                continue;
            }else{
                break;
            }
        }
             
        if(lastIdx == 0){
            lastIdx = lastIdx;
        }
        
        final double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            CspProcess.bInTransition = false;
        }

        
        final int totalCons = 8;
        double t[] = new double[totalCons];

        t[0] = x[0] + x[1]*(x[0] + x[2]) + x[3]*(x[2] + x[4]) + x[5]*(x[4] + x[6]) - x[7]*((1/8) - x[6]);// = 0; 
        t[1] = x[1] + x[2]*(x[0] + x[4]) + x[3]*(x[1] + x[5]) + x[4]*x[6] - x[7]*((2/8) - x[5]);//=0
        t[2] = x[2]*(1 + x[5]) + x[3]*(x[0] + x[6]) + x[1]*x[4] - x[7]*((3/8) - x[4]);//=0
        t[3] = x[3] + x[0]*x[4] + x[1]*x[5] + x[2]*x[6] - x[7]*((4/8) - x[3]);//=0
        t[4] = x[4] + x[0]*x[5] + x[1]*x[6] - x[7]*((5/8) - x[2]);//=0
        t[5] = x[5] + x[0]*x[6] - x[7]*((6/8) - x[1]);//=0
        t[6] = x[6] - x[7]*((7/8) - x[0]);//=0
        t[7] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] +1; //= 0;
                
        for (int i = 0; i < t.length; i++) {
            consFn.add(new Constraint(t[i]));            
        }
                
        for (int i = 0; i < consFn.size(); i++) {
            if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
            }else{
                consFn.get(i).isViolated = true;
                violation.add(i);
            } 
        }
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = t.length; i < t.length+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }
        
    }
    
    private static void hs109(final int tm) throws SolutionFoundException{
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 9
        // total constraints = 11+9
		
        // 0<= x[i] <= inf; i = 0,1 
        // -.55 <= x[i] <= .55; i = 2,3
	// 196 <= x[i] <= 252; i = 4, 5, 6,
        // -400 <= x[i] <= 800; i = 7, 8
        
        curFitness = 0.0;
        violation.clear();
        
        final int minPow = -1;
        final int maxPow = 20;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow-minPow)+1;
        
//        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
//        for (int i = 1; i < total_divisions; i++) {
//            if(lastIdx>grpLenght*i){
//                grp++;
//                continue;
//            }else{
//                break;
//            }
//        }
             
        grp = CspProcess.negFeasibleRange;//Math.min(grp,CspProcess.negFeasibleRange);         
        
        double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            dt = MIN_DT;
            CspProcess.bInTransition = false;
        }
        
        
        
        
        double a = 50.176;
        double b1 = 0.25;

        double b = Math.sin(b1);
        double c = Math.cos(b1);

        double best_val_found = 5327.541669;
        double eps = 532.7541669; //10% error	

        final int tabuSize = tabuConstraints.size();
        final int totalFunctionalConstraints = 11;     
        double t[] = new double[totalFunctionalConstraints+tabuSize];
//        double t[] = new double[1];

        t[0] = best_val_found + eps - (3.0 * x[0] + Math.pow(10.0,-6) * Math.pow(x[0],3) + 2 * x[1] + 
                0.522074*Math.pow(10.0,-6) * Math.pow(x[1],3)); //>=0
        t[1] = x[3] - x[2] + 0.55;// >= 0;
        t[2] = x[2] - x[3] + 0.55;// >= 0;
        t[3] = 2250000.0 - Math.pow(x[0],2) - Math.pow(x[7],2);// >= 0;
        t[4] = 2250000.0 - Math.pow(x[1],2) - Math.pow(x[8],2);// >= 0;
        
        t[5] = x[4] * x[5] * Math.sin(-x[2] - 0.25) + x[4] * x[6] * Math.sin(-x[3] - 0.25) +
                2 * b * Math.pow(x[4],2) - a * x[0] + 400.0 * a;// = 0;
        t[6] = x[4] * x[5] * Math.sin(x[2] - 0.25) + x[5] * x[6] * Math.sin(x[2] - x[3] - 0.25) +
                2 * b * Math.pow(x[5],2) - a * x[1] + 400.0 * a;// = 0;
        t[7] = x[4] * x[6] * Math.sin(x[3] - .25) + x[5] * x[6] * Math.sin(x[3] - x[2] - .25) +
                2 * b * Math.pow(x[6],2) + 881.779 * a;// = 0;
        t[8] = a * x[7] + x[4] * x[5] * Math.cos(-x[2] - 0.25) +
                x[4] * x[6] * Math.cos(-x[3] - 0.25) - 200.0 * a - 2 * c * Math.pow(x[4],2) +
                0.7533*Math.pow(10.0,-3.0) * a * Math.pow(x[4],2);// = 0;
        t[9] = a * x[8] + x[4] * x[5] * Math.cos(x[2] - .25) +
                x[5] * x[6] * Math.cos(x[2] - x[3] - 0.25) - 2.0 * c * Math.pow(x[5],2) +
                0.7533*Math.pow(10.0,-3.0) * a * Math.pow(x[5],2) - 200.0 * a; // = 0;
        t[10] = x[4] * x[6] * Math.cos(x[3] - 0.25) + x[5] * x[6] * Math.cos(x[3] - x[2] - 0.25) -
                2 * c * Math.pow(x[6],2) + 22.938 * a + 
                0.7533*Math.pow(10.0,-3.0) * a *Math.pow(x[6],2); // = 0;
        
        double sum=0.0;
        if(tabuSize>0){
            for (int cons = 0; cons < tabuConstraints.size(); cons++) {
                ArrayList<Double> centers = tabuConstraints.get(cons);
                if(x.length != centers.size()-1){ //last one is radius
                    throw new SolutionFoundException("NO SOLUTION! BUG!!!! Center and varX size must be same. Check your code");
                }
                sum = 0.0;
                //This is equation of a hyper-sphere
                for (int i = 0; i < x.length; i++) {
                    sum += Math.pow(x[i]-centers.get(i),2);                                        
                }
                //hyper-sphere equation. (x1-c1)^2 + (x2-c2)^2 + .... +(xn-cn)^2 - r^2 >= 0
                t[totalFunctionalConstraints+cons] = sum - Math.pow(centers.get(centers.size()-1),2); //>=0
            }
        }        
        
        //Short term memory - for more than 5 tabus start deleting/forgetting the oldest one.
        final int SHORT_TERM_MEMORY = 5;
        for (int i = SHORT_TERM_MEMORY; i < tabuSize; i++) {
            t[totalFunctionalConstraints+(i-SHORT_TERM_MEMORY)] = 1.0; //to make it >=0
        }
        
        final int USED_CONS = Math.min(tm,t.length);

        int consOrder[] = new int[100];//{0,1,2,3,4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, };//{0,3,4,5,6,9,10,8,1,2,7};
        Arrays.fill(consOrder, -1);
        int consOrderFncons[] =  {7,2,1,8,10,9,6,0,3,4,5};//{5,4,3,0,6,9,10,8,1,2,7}; //
        for (int i = 0; i < consOrderFncons.length; i++) {
            consOrder[i] = consOrderFncons[i];
            
        }        
        for (int i = totalFunctionalConstraints; i < t.length; i++) {
            consOrder[i] = i;            
        }
//        int consOrder[] = {0};
        
        
        for (int i = 0; i < USED_CONS; i++) {
            consFn.add(new Constraint(t[consOrder[i]]));            
        }
        
  
        int idx; 
        for (int i = 0; i < USED_CONS; i++) {
            idx = consOrder[i];
            
            if((idx>=0 && idx <=4) || idx >= totalFunctionalConstraints){
                if(consFn.get(i).val >= 0.0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                } 
            }else{
                if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                } 
            }            
        }
          
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = USED_CONS; i < USED_CONS+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-USED_CONS) || 
            consFn.get(i).val>userInput_.maxVals.get(i-USED_CONS)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }		
    }
    
    
//    private static void h77(){
//        Double x[] = new Double[userInput_.totalDecisionVars];
//        vals.toArray(x);
//        // total vars (n) = 5
//        // total constraints = 3+5
//        // 0 <= xi <= 4, i = 0..4
//
//        
//        final double best_val_found = 0.2415051288;
//        final double eps = 0.002415051288;
//        final double s2 = Math.sqrt(2.0);
//        
//        curFitness = 0.0; //It is CSP problem
//        violation.clear();
//        
//        
//        final int minPow = -3;
//        final int maxPow = 2;
//        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
//        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
//        final int total_divisions = Math.abs(maxPow-minPow)+1;
//        
//        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
//        int grp = 0;
//        for (int i = 1; i < total_divisions; i++) {
//            if(lastIdx>grpLenght*i){
//                grp++;
//                continue;
//            }else{
//                break;
//            }
//        }
//               
//        final double dt = MAX_DT*Math.pow(10, -grp); 
//        
//        if(dt > MIN_DT){
//            CspProcess.bInTransition = true;
//        }else{
//            CspProcess.bInTransition = false;
//        }
//        
//        
//        double t[] = new double[3];
//        t[0] = best_val_found + eps - (2*x[0]*(x[0]-x[1]-1) + 1 + x[1]*x[1] + 
//                (x[2] - 1)*(x[2] - 1) + Math.pow(x[3] - 1.0,4.0) + Math.pow(x[4] - 1.0,6.0));//>=0;
//        t[1] = (x[0]*x[0]) * x[3] + Math.sin(x[3] - x[4]) - 2.0 * s2;//=0;
//        t[2] = x[1] + Math.pow(x[2]*1.0,4.0)*(x[3]*x[3]) - 8 - s2;//=0;
//        
//        for (int i = 0; i < t.length; i++) {
//            consFn.add(new Constraint(t[i]));            
//        }  
//        
//        
//        
//        if(consFn.get(0).val >= 0){
//            CSPsols.get(lastIdx).set(0, (ArrayList<Double>)vals.clone());   
//        }else{
//            consFn.get(0).isViolated = true;
//            violation.add(0);
//        }
//        
//        for(int i = 1; i< t.length; i++){
//            if(dt - Math.abs(consFn.get(i).val) >= 0.0){
//                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
//            }else{
//                consFn.get(i).isViolated = true;
//                violation.add(i);
//            } 
//        }
//    
//        for (int i = 0; i < x.length; i++) {
//            consFn.add(new Constraint(x[i])); 
//        }
//                
//        for (int i = t.length; i < t.length+x.length; i++) {            
//            if(consFn.get(i).val <userInput_.minVals.get(i-t.length) || 
//                    consFn.get(i).val>userInput_.maxVals.get(i-t.length)){
//                consFn.get(i).isViolated = true;
//                violation.add(i);
//            }else{
////                 if(CSPsols.size() != maxCSPsolsSize-1 && CSPsols.get(lastIdx).get(i).isEmpty())
//                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
//            }
//        }
//        
//    }
    private static void h77(final int tm) throws SolutionFoundException{
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        // total vars (n) = 5
        // total constraints = 3+5
        // 0 <= xi <= 4, ii = 0..4

        
        final double best_val_found = 0.2415051288;
        final double eps = 0.002415051288;
        final double s2 = Math.sqrt(2.0);
        
        curFitness = 0.0; //It is CSP problem
        violation.clear();
        
        
        final int minPow = -3;
        final int maxPow = -1;
        final double MIN_DT = Math.pow(10, minPow); //1.0E-1;
        final double MAX_DT = Math.pow(10, maxPow); //1.0E2;
        final int total_divisions = Math.abs(maxPow -minPow)+1;
        
//        int grpLenght = maxCSPsolsSize/total_divisions;//need floor int
        int grp = 0;
//        for (int i = 1; i < total_divisions; i++) {
//            if(lastIdx>grpLenght*i){
//                grp++;
//                continue;
//            }else{
//                break;
//            }
//        }
        
        grp = CspProcess.negFeasibleRange;//Math.min(grp,CspProcess.negFeasibleRange);
        
               
        double dt = MAX_DT*Math.pow(10, -grp); 
        
        if(dt > MIN_DT){
            CspProcess.bInTransition = true;
        }else{
            dt = MIN_DT;
            CspProcess.bInTransition = false;
        }
        
        final int tabuSize = tabuConstraints.size();
        final int totalFunctionalConstraints = 3;     
        double t[] = new double[totalFunctionalConstraints+tabuSize];
//        double t[] = new double[1];
        

        t[0] = best_val_found + eps - (2*x[0]*(x[0]-x[1]-1) + 1 + x[1]*x[1] + 
                (x[2] - 1)*(x[2] - 1) + Math.pow(x[3] - 1.0,4.0) + Math.pow(x[4] - 1.0,6.0));//>=0;
        t[1] = (x[0]*x[0]) * x[3] + Math.sin(x[3] - x[4]) - 2.0 * s2;//=0;
        t[2] = x[1] + Math.pow(x[2]*1.0,4.0)*(x[3]*x[3]) - 8 - s2;//=0;
        
        double sum=0.0;
        if(tabuSize>0){
            for (int cons = 0; cons < tabuConstraints.size(); cons++) {
                ArrayList<Double> centers = tabuConstraints.get(cons);
                if(x.length != centers.size()-1){ //last one is radius
                    throw new SolutionFoundException("NO SOLUTION! BUG!!!! Center and varX size must be same. Check your code");
                }
                sum = 0.0;
                //This is equation of a hyper-sphere
                for (int i = 0; i < x.length; i++) {
                    sum += Math.pow(x[i]-centers.get(i),2);                                        
                }
                //hyper-sphere equation. (x1-c1)^2 + (x2-c2)^2 + .... +(xn-cn)^2 - r^2 >= 0
                t[totalFunctionalConstraints+cons] = sum - Math.pow(centers.get(centers.size()-1),2); //>=0
            }
        }
        
        
        final int USED_CONS = Math.min(tm,t.length);
      
        int consOrder[] = new int[100];//{0,1,2,3,4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, };//{0,3,4,5,6,9,10,8,1,2,7};
        Arrays.fill(consOrder, -1);
        int consOrderFncons[] = {0,1,2};
        for (int i = 0; i < consOrderFncons.length; i++) {
            consOrder[i] = consOrderFncons[i];
            
        }        
        for (int i = totalFunctionalConstraints; i < t.length; i++) {
            consOrder[i] = i;            
        }

//        int consOrder[] = {0};
        
        for (int i = 0; i < USED_CONS; i++) {
            consFn.add(new Constraint(t[consOrder[i]]));            
        } 
        
        int idx; 
        for (int i = 0; i < USED_CONS; i++) {
            idx = consOrder[i];
            
            if(idx==0 || idx >= totalFunctionalConstraints){       
                if(consFn.get(i).val >= 0){
                    CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                }
            }else{
                if(dt - Math.abs(consFn.get(i).val) >= 0.0){
                        CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());   
                }else{
                    consFn.get(i).isViolated = true;
                    violation.add(i);
                } 
            }
        }
    
        
        
        for (int i = 0; i < x.length; i++) {
            consFn.add(new Constraint(x[i])); 
        }
                
        for (int i = USED_CONS; i < USED_CONS+x.length; i++) {            
            if(consFn.get(i).val <userInput_.minVals.get(i-USED_CONS) || 
            consFn.get(i).val>userInput_.maxVals.get(i-USED_CONS)){
                consFn.get(i).isViolated = true;
                violation.add(i);
            }else{
                CSPsols.get(lastIdx).set(i, (ArrayList<Double>)vals.clone());  
            }
        }        
    }
    
    private static void objFunction_h74mod(){
        
        fitness_.clear();
        

        violation.clear();
        Double x[] = new Double[userInput_.totalDecisionVars];
        vals.toArray(x);
        
        double a = 0.55;  
        double c1 = Math.pow(2.0/3.0,-6);

//	cons1 {i in 1..2}: x[i] >= 0;
//cons2 {i in 1..2}: x[i] <= 1200;
//cons3 {i in 3..4}: x[i] >= -a;
//cons4 {i in 3..4}: x[i] <= a;


        if(x[0] * (3 + Math.pow(1.0,-6.0)*x[0]*x[0]) + x[1] * (2 + Math.pow(c1*x[1],2.0)) > 5126.49811 + 51.2649811);
            violation.add(1);

        if(x[3] - x[2] + a <0)
            violation.add(2);

        if(x[2] - x[3] + a <0)
            violation.add(3);

        if(1000 * Math.sin(-x[2] - 0.25) + 1000 * Math.sin(-x[3]-0.25) + 894.8 - x[0] != 0)
            violation.add(4);	

        if(1000 * Math.sin(x[2] - 0.25) + 1000 * Math.sin(x[2]-x[3]-0.25) + 894.8 - x[1] != 0)
            violation.add(5);

        if(1000 * Math.sin(x[3] - 0.25) + 1000 * Math.sin(x[3]-x[2]-0.25) + 1294.8 != 0)
            violation.add(6);
        
    
//    for (Double cf : consFn) {
//        if(cf<0){
//            position.add(-1);
//        }else if (cf>0){
//            position.add(1);
//        }else{
//            position.add(0);
//        }            
//    }
    
        fitness_.clear();
        if(maxCSPval>= Double.MAX_VALUE)
            fitness_.add(violation.size()*1.0);
        else
            fitness_.add(0+violation.size()*maxCSPval);
    }
    
    
        private static void objFunction_deb(){
    //<< deb problem

        ArrayList<Double> x = vals;

        violation.clear();

        //if(x.get(1) + 5.3*x.get(0) < 6){
        if(x.get(1) + 9*x.get(0) < 6){
            violation.add(0);
        }else{
            CSPsols.get(lastIdx).set(0, (ArrayList<Double>)vals.clone());
        }

        //if(-x.get(1) + 1.8*x.get(0) < 1){
        if(-x.get(1) + 9*x.get(0) < 1){
            violation.add(1);
        }else{
            CSPsols.get(lastIdx).set(1, (ArrayList<Double>)vals.clone());
        }

        if(x.get(0) > 1){
            violation.add(2);
        }else{
            CSPsols.get(lastIdx).set(2, (ArrayList<Double>)vals.clone());
        }
        
        
        fitness_.clear();
        if(maxCSPval>= Double.MAX_VALUE)
            fitness_.add(violation.size()*1.0);
        else
            fitness_.add(0+violation.size()*maxCSPval);
    //>> deb problem
    }
    
    
}
