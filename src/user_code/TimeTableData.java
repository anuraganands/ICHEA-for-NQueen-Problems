///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package user_code;
//
//import csp.Chromosome;
//import csp.Element;
//import csp.ExternalData;
//import csp.Idx2D;
//import csp.MyException;
//import csp.MyRandom;
//import csp.SolutionFoundException;
//import java.io.File;
//import java.io.FileInputStream;
//import java.io.FileNotFoundException;
//import java.io.ObjectInputStream;
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.HashSet;
//import java.util.Scanner;
//import java.util.StringTokenizer;
//import org.jdesktop.application.Application;
//
//
///**
// *
// * @author s425770
// */
//// It should be singleton object
//    
//public class TimeTableData extends ExternalData{
//    private ArrayList<ArrayList<Integer>> courses_; //10 courses_ 0-9
//    private int conflictMatrix_[][];
//    private ArrayList<Element> conflictOrder; //it contains sorted order of courses from most constrained to least constrained.
//    int coeffConflictVsCourse[]; 
//    private int maxSlots;//maximum grouped constraint satisfaction allowed
//    private boolean bReadAllFromFile;
//    int STUDENT_SIZE;
//    int COURSE_SIZE;
//    boolean bdebug;
//
//    
//    public TimeTableData(String fileName, int populaion, int generation, int curPref, 
//            int prevPref, boolean saveChromes, int solutionBy, Class t) 
//            throws InstantiationException, IllegalAccessException, MyException {
//        super(fileName, populaion, generation, curPref, prevPref, saveChromes, solutionBy, t);        
//        maxSlots = 0;
//        initializeCounter = 0;
//        this.courses_ = new ArrayList<ArrayList<Integer>>();
//        
//        bReadAllFromFile = false;
//        conflictOrder = new ArrayList<Element>();
//        bdebug = false;
//        
//        this.readData();
//
//    }
//
//    @Override
//    public void readData() throws MyException {
//        //Students ID starting from 0 to max val
//        //Courses ID starting frmo 0 to max val
//
//        ArrayList<Integer> arrayList;
//        
//        ArrayList<Integer> students1;
//        ArrayList<Integer> students2;
//        int conflictStudent;
//                
//        Scanner dataFile;
//        try {
//            dataFile = new Scanner(new File(fileName_));
//        
//        
//            while (dataFile.hasNext("#")){                
//                System.out.println("got #####");
//                dataFile.nextLine(); //ignore comments  
//            }
//
//            STUDENT_SIZE = dataFile.nextInt();
//            dataFile.nextLine(); //ignore comments
//
//            COURSE_SIZE = dataFile.nextInt();
//            dataFile.nextLine(); //ignore
//
//            maxSlots = dataFile.nextInt();
//            dataFile.nextLine(); //ignore             
//
//            for (int i = 0; i < COURSE_SIZE; i++) {
//                arrayList = new ArrayList<Integer>();
//                courses_.add(arrayList);  
//            }
//
//            StringTokenizer str;
//            Integer tempCourse;
//            int tempStudentID = 0;
//            int sz;
//            String tempStr;
//
//            while(dataFile.hasNext()){    
//                tempStr = dataFile.nextLine();
//                str = new StringTokenizer(tempStr," ");
//                sz = str.countTokens();
//                for (int i = 0; i < sz; i++){
//
//                    tempStr = str.nextElement().toString();                
//                    tempCourse = Integer.valueOf(tempStr);
//                    courses_.get(tempCourse-1).add(tempStudentID);
//                }  
//                if(sz>0)
//                    tempStudentID++;                
//            }
//        
//            //userInput_.fileData = true;
//            userInput_.totalConstraints = courses_.size();
//            userInput_.totalDecisionVars = userInput_.totalConstraints;
//            userInput_.totalObjectives = 1;
//            userInput_.bHasConstraintPreferences = true;
//            userInput_.bWeighted = false; 
//            for (int i = 0; i < userInput_.totalDecisionVars; i++) {
//                userInput_.minVals.add(0.0);
//                userInput_.maxVals.add(userInput_.totalConstraints*1.0-1);            
//            }  
//
//            userInput_.validateData();
//            userInput_.domainVals = new  ArrayList<ArrayList<Double>>();
//            for (int i = 0; i < userInput_.totalDecisionVars; i++) {
//                userInput_.domainVals.add(new ArrayList<Double>());
//                for (Double j = userInput_.minVals.get(i); j <= userInput_.maxVals.get(i); j++) {
//                    userInput_.domainVals.get(i).add(j);
//                }
//            }        
//            userInput_.doMutation = true;
//            immunitySize = 0;//userInput_.totalDecisionVars; //this immunity is different from immunity used in chromosome.
//        
//            conflictMatrix_ = new int[COURSE_SIZE][COURSE_SIZE];
//            coeffConflictVsCourse = new int[COURSE_SIZE];
//            for (int i = 0; i < courses_.size(); i++) {
//                students1 = courses_.get(i);
//                for (int j = i; j < courses_.size(); j++) {
//                    if(i == j){
//                        conflictMatrix_[i][j] = students1.size();
//                    }
//                    students2 = courses_.get(j);     
//                    conflictStudent = 0;
//                    for (Integer astudent : students1) {
//                        if(students2.contains(astudent)){
//                            conflictStudent++;
//                        }
//                    }
//                    conflictMatrix_[i][j] = conflictStudent;
//                    conflictMatrix_[j][i] = conflictStudent;
//                }
//            }
//            
//            int total;
//            for (int i = 0; i < COURSE_SIZE; i++) {
//                total = 0;
//                for (int j = 0; j < COURSE_SIZE; j++) {
//                    if(conflictMatrix_[i][j]>0)
//                        total += 1;                    
//                }
//                conflictOrder.add(new Element(total*1.0, i, Element.DESCENDING));
//            }
//            
//            Collections.sort(conflictOrder);
//            
//            for (int wt = 0; wt < conflictOrder.size(); wt++) {
//                coeffConflictVsCourse[conflictOrder.get(wt).idx] = conflictOrder.size() - wt;                
//            }    
//
//        } catch (FileNotFoundException fnfe) {
//            System.err.println("Data file not found.");
//            Application.getInstance().exit();
//        }
//        //debug();
//    }
//
////    @Override
////    public ArrayList<Double> negateVal(ArrayList<Double> vals) {
////        throw new UnsupportedOperationException("Not supported yet.");
////    }
//
//    
//    /**
//     * Don't try to use fitness value here becuase it is corrupted...
//     * @param chromeConstraints
//     * @param bShowProgress
//     * @return 
//     */
//    @Override
//    public boolean getForcedCSPsol(ArrayList<ArrayList> chromeConstraints, boolean bShowProgress) {
//            ArrayList<Double> vals = new ArrayList<Double>();
//            ArrayList<Double> fitness = new ArrayList<Double>(); //NOTE fitness will be very corrupted.... 
//            double tmpVal;
//            int tmpIdx;
//            
//            chromeConstraints.clear();
//            ArrayList<Double> noGood = new ArrayList<Double>();
//            Idx2D[] valVsConstIdx_ = new Idx2D[userInput_.totalConstraints];
//            
//            for (int wt = 0; wt < maxSlots; wt++) {
//                chromeConstraints.add(new ArrayList<Integer>());            
//            }
//            
//            int bs, as;
//            Double badval;
//            Integer temp;
//            int valSize;
//            int ConstSize;
//            
//            ArrayList<Integer> idxList = MyRandom.randperm(0,conflictOrder.size()-1);
//            
//            for (int i = 0; i < conflictOrder.size(); i++) {
//            //for(int i = 0; i<idxList.size();i++){
//                tmpIdx = i;//idxList.get(i);
//                if(!fitness.isEmpty()){
//                    fitness.add(Integer.MAX_VALUE*1.0);
//                    fitness.set(0, Integer.MAX_VALUE*1.0);// fitness is already corrupted....
//                }
//                tmpVal = conflictOrder.get(tmpIdx).idx*1.0;
//                vals.add(tmpVal);
//                if(vals.subList(0, vals.size()-1).contains(tmpVal)){
//                    vals.remove(vals.size() - 1);
//                    continue;
//                }
//                
//                bs = vals.size();
//                ObjectiveFnAppend(vals, fitness, chromeConstraints, noGood, valVsConstIdx_);
//                as = vals.size();
//                
//                if(as < bs){//not added
//                    badval = tmpVal;
//                    i = -1; //start again... but these will be ignored by append fn
//                    
//                    int j = MyRandom.randperm(0, chromeConstraints.size()-1).get(0);
//                    for (int k = 0; k < chromeConstraints.get(j).size(); k++) {
//                        //System.out.println(chromeConstraints.get(j).get(k));
//                        temp = (Integer)chromeConstraints.get(j).get(k);
//                        if(isViolated(temp, badval.intValue())){
//                            chromeConstraints.get(j).remove(k);
//                            vals.remove(temp*1.0);
//                            k--;
//                        }
//                    } 
//                    vals.add(badval);
//                    ObjectiveFnAppend(vals, fitness, chromeConstraints, noGood, valVsConstIdx_);
//                } 
//                
//                if(i%5 == 0 && bShowProgress){
//                    System.out.println(chromeConstraints+"\n\n");
//                }
//            }
//            valSize = vals.size();
//            ConstSize = 0;
//            for (ArrayList grp : chromeConstraints) {
//                for (Object obj : grp) {
//                    ConstSize++;
//                }
//            }
//            fitness.clear(); //note fitness val is corrupted...
//            if(ConstSize == valSize && valSize == userInput_.totalConstraints){
//                return true;
//            }else
//                return false;
//    }
//
//        /**
//     * Don't try to use fitness value here becuase it is corrupted...
//     * @param chromeConstraints
//     * @param bShowProgress
//     * @return 
//     */
//    @Override
//    public void tryForcedCSPsolUpdate(final ArrayList<Double> vals, final ArrayList<Double> fitness, 
//    final ArrayList<ArrayList> chromeConstraints, final ArrayList<Double> noGood,
//    final Idx2D[] valVsConstIdx_,  boolean bShowProgress){
//        //ArrayList<Double> vals = new ArrayList<Double>();
//        //ArrayList<Double> fitness = new ArrayList<Double>(); //NOTE fitness will be very corrupted.... 
//        double tmpVal;
//        int tmpIdx;
//
//        ArrayList<Double> localNoGoods = new ArrayList<Double>();
//        double randNoGood;
//        int randListIdx;
//        ArrayList<Integer> randList;
//
//        
//        for (int i = 0; i < userInput_.totalConstraints; i++) {
//            localNoGoods.add(i*1.0);
//        }
//        
//        for (int i = 0; i < vals.size(); i++) {   
//            localNoGoods.remove(vals.get(i));        
//        }
//        
//        if(localNoGoods.isEmpty()){
//            return;
//        }
//
//        randNoGood = localNoGoods.get(MyRandom.randperm(0, localNoGoods.size()-1).get(0));
//
//        //It is assumed that all groups are conflicted with localNoGoods... duh... that's the definition of localNoGoods...
//        randListIdx = MyRandom.randperm(0, maxSlots-1).get(0);
//        randList = chromeConstraints.get(randListIdx); //reference to list.
//
//
//        for (int j = 0; j < randList.size(); j++) {
//            if(isViolated(randList.get(j).intValue(), (int)randNoGood)){
//                valVsConstIdx_[randList.get(j)] = new Idx2D();
//                randList.remove(j);
//                j--;                
//            }
//        }
//        randList.add((int)randNoGood);
//        
//        chromeConstraints.set(randListIdx, randList);
//        
//        ObjectiveFnRefresh(chromeConstraints, fitness, vals, valVsConstIdx_);
//        
////        for (int i = 0; i < randList.size(); i++) {
////            valVsConstIdx_[randList.get(i)].position = i;            
////        }
//        
//        //ObjectiveFnAppend(vals, fitness, chromeConstraints, noGood, valVsConstIdx_);
//           
//    }
//    
//    
//    
//    @Override
//    public ArrayList<Chromosome> initializeExternalChrmosomes(final int population) throws SolutionFoundException{ //, final boolean bInitialStage) {
//        Chromosome chrome;
//        ArrayList<Chromosome> newChromes = new ArrayList<Chromosome>();
//        int modVal;
//        //int prevPref;
//        
//        if(userInput_ == null)
//            throw new UnsupportedOperationException("User input not initialized");
//    
//        if(!bReadAllFromFile && this.getPrevPref() >= 0){
//            try {           
//                FileInputStream fis = new FileInputStream("partial_solutions_pref_" + this.getPrevPref() + ".ichea");
//                ObjectInputStream ois = new ObjectInputStream(fis);
//                newChromes = (ArrayList<Chromosome>)ois.readObject();
//                
//                this.immunitySize = newChromes.get(0).getValsCopy().size(); //I hope only satisfaction is used...
//                this.pStart = newChromes.get(0).getExternalData().getpStart();
//                
//                for (Chromosome nc : newChromes) {
//                    nc.setExternalData(this);
//                }
//                ois.close();
//            }
//            catch(Exception e) {
//                System.out.println("Exception during deserialization: " + e);
//                System.exit(0);
//            }
//
//            bReadAllFromFile = true; // rest are discarded ... as in next call to initializeExternalChrmosomes data will not be read from a file
//            if(newChromes.size()>population){
//                newChromes = (ArrayList<Chromosome>)newChromes.subList(0, population);
//                return newChromes;
//            }
//        }
//        
////        if(bInitialStage){
//            modVal = userInput_.totalConstraints;
////        }else{
////            modVal = maxSlots;
////        }
//        
//        for (int i = 0; newChromes.size() < population; i++) { 
//            chrome = new Chromosome(userInput_.solutionBy, this); 
//            for (int j = 0; j < maxSlots; j++) {                
//                chrome.appendVal(conflictOrder.get(j).idx, Double.MAX_VALUE);
//            }
//   
//            chrome.appendVal(conflictOrder.get(initializeCounter++ % modVal).idx, Double.MAX_VALUE );
//            newChromes.add(chrome);                     
//        }
//                  
//        if(bdebug){
//            userInput_.bWeighted = true;
////            debug(newChromes.get(0).getSatisfaction(), newChromes.get(0).fitness_);  
////            System.out.println("The input solution is: " + newChromes.get(0));
////            System.out.println("The soulution is: "+newChromes.get(0).fitness_.get(0)/STUDENT_SIZE);
//            Application.getInstance().exit();
//        }
//
//        return newChromes;
//    }
//
//    /**
//     * Check violation in 2 given courses_. It checks if there is overlap of 2
//     * or more students.
//     * @param obj1 - First course index - note index starts from 0
//     * @param obj2 - Second course index - note index starts from 0
//     * @return if there is overlap of 2 or more students that means there is
//     * a violation and true is returned, otherwise false is returned.
//     */
//
//    
//    private int penaltyFn(int dist){
//        int weight = Integer.MAX_VALUE;
//        if(dist == 0){
//            weight = 32;//Integer.MAX_VALUE; //16;
//        }else if (dist == 1){
//            weight = 16; //8;
//        }else if (dist == 2){
//            weight = 8;//4;
//        }else if (dist == 3){
//            weight = 4;//2;
//        }else if (dist == 4){
//            weight = 2; //1;
//        }else if(dist > 4){
//            weight = 1;//0;
//        }else{
//            throw new UnsupportedOperationException("dist must be > 0.");
//        }
//        return weight;
//    }
//    
//    /**
//     * 
//     * @param dist maximum value of dist = levels - 1;
//     * @return 
//     */
//    private int prefFn(int dist){
//        int levels = 5;
//        int pref = Integer.MAX_VALUE;
//        
//        if(dist > 0){
//            pref = levels - dist;
//            if(pref < 0)
//                pref = 0;
//        }else{
//            ;
//            //throw new UnsupportedOperationException("dist cannot be more than max level defined.");
//        }
//        
//        return pref;
//    }
//
//    /**
//     * 
//     * @param fnVal it is prefFn or penaltyFn
//     * @param satisfiedCons - total satisifactions
//     * @return fitness value
//     */
////    private double curOjbectiveFn(final int maxPref, final int prefVal, final int satisfiedCons){
////        //prefVal = (prefVal+1)*(userInput_.totalConstraints - satisfiedCons+1)-1;
////        //return prefVal*1.0;
////        /// f(Lenght, pref) = f(k,p) = (2.MaxP+1)^(L-l).(p+1) 
////        double val = Math.pow(2.0*maxPref+1,userInput_.totalConstraints*1.0-satisfiedCons)*(prefVal + 1.0);
////        //return val/Math.pow(2.0*maxPref, userInput_.totalConstraints/2.0);
////        return val;
////    }
//    
//
//    /**
//     * 
//     * @param p current preference = penalty*effected_students
//     * @param l total constraints satisfied so far (or total courses so far)
//     * @param fitness
//     * @return 
//     */
//    @Override
//    public double fitnessValWeightBased(final int penalty, ArrayList<Double> fitness){        
//        //double ifitness = curStudentSize*penaltyFn(dist); 
//        final int L = userInput_.totalConstraints;
//        final double totalWorstCaseVal;
//        final double individualWorstCaseVal;
//        final double maxPenalty = penaltyFn(0);
//        
//                
//        totalWorstCaseVal = 1.0*L*(L-1)*STUDENT_SIZE*maxPenalty;
//        individualWorstCaseVal = 1.0*STUDENT_SIZE*(L-1)*maxPenalty;
//        
//        double curFitness;
//        
//        if(fitness.isEmpty()){
//            fitness.add(totalWorstCaseVal);                 
//        }  
//
//        curFitness = fitness.get(0)-individualWorstCaseVal;
//        curFitness+= penalty;    
//        
//        if(fitness.get(0)<0 || curFitness < 0){
//            System.out.println("eee kaisey....");
//        }  
//
//        fitness.set(0,curFitness); //one just added          
//        
//
//        
//        return fitness.get(0);        
//    }
//    private double curOjbectiveFn(final int D, final int p, final int constrainedWt, ArrayList<Double> fitness){
//        if(userInput_.bWeighted){
//            return fitnessValWeightBased(p, fitness);
//        }else{
//            return fitnessValPrefBased(D, p, constrainedWt, fitness);
//        }
//    }
//
//    @Override
//    protected int maxPref() { 
//        int D = prefFn(0);
//        if(userInput_.bWeighted){
//            //return penaltyFn(0);
//            //return 2*D*STUDENT_SIZE*penaltyFn(0);
//            return Integer.MAX_VALUE;   //??? it should be 2*D*Student_SIze*penalty(0)???         
//        }else{
//            return 5;//prefFn(0);
//        }
//    }
//
//    @Override
//    protected void ObjectiveFnRefresh(ArrayList<ArrayList> constraints, ArrayList<Double> fitness_, 
//    ArrayList<Double> vals, Idx2D[] valVsConstIdx_) {
//        debug(constraints, fitness_);
//        
//        vals.clear();
//        for (ArrayList grp : constraints) {
//            for (Object obj : grp) {
//                vals.add(1.0*((Integer)obj));
//            }
//        }
//        
//        Idx2D idx2D;
//        int col;
//        int pos;
//        //int i;
//         
//        for (int j = 0; j < valVsConstIdx_.length; j++) {
//            valVsConstIdx_[j] = new Idx2D();
//        }
//
//        col = -1;
//        //i = 0;
//        for (ArrayList<Integer> list : constraints) {
//            col++;
//            pos = -1;
//            for (Integer val : list) {
//                pos++;
////                idx2D = new Idx2D();
////                idx2D.col = col;
////                idx2D.position = pos;
////                valVsConstIdx_[val] = idx2D;
//                valVsConstIdx_[val].col = col;
//                valVsConstIdx_[val].position = pos;
//               // i++;
//            }
//        }
//
//    }
//
//    
//    
//    
//    /**
//     * 
//     * @param vals - NOTE vals must be unique.
//     * @param fitness - fitness[0] = fitness val, fitness[1] = Sum of preferences
//     * fitness[2] = worst preference encountered so far.
//     * @param chromeConstraints 
//     */
//    @Override
//    protected void ObjectiveFnAppend(final ArrayList<Double> vals, final ArrayList<Double> fitness, 
//    final ArrayList<ArrayList> chromeConstraints, final ArrayList<Double> noGood,
//    final Idx2D[] valVsConstIdx_) {
//        //All values will be accepted..
//        ArrayList<Integer> tempAL;
//        ArrayList<Integer> Prefs = new ArrayList<Integer>();
//        ArrayList<Integer> conflictStudents = new ArrayList<Integer>();
//        ArrayList<Integer> possibleLoc = new ArrayList<Integer>();
//        final int course = vals.get(vals.size()-1).intValue();
//        final int MAX_DIST = 5;
//        int MAX_PREF_WT;
//        int MIN_PREF_WT;       
//        int totalConflictStudents;
//        int worstPref = 0;
//        boolean bPlaceFound = false;
//        int tempDist;  
//        int hardConstraintViolated = 0;
//        Idx2D idx2D;
//        int minIdx = -1;
//
//
//        MIN_PREF_WT = 0;
//        MAX_PREF_WT = 5; //prefFn(0);
//
//
//        if(vals.size() == 0)
//            return;
//        
//        
//       //-wt think eela nikalek pari... 
//        //expensive bhi hai aur ... setfn to deal kare hai... lekin appendfn nahi.... 
//        
//        if(vals.subList(0, vals.size()-1).contains(course*1.0)){
//            vals.remove(vals.size() - 1);//last element
//            //fitness and satisfactions not changed.
//            return;
//        }
//        
//        if(vals.size() == 1){
//            chromeConstraints.clear();
//            for (int i = 0; i < maxSlots; i++) {
//                chromeConstraints.add(new ArrayList<Integer>());                
//            }
//
//            chromeConstraints.get(0).add(vals.get(0).intValue());
//            
//            idx2D = new Idx2D();
//            idx2D.col = 0;
//            idx2D.position = 0;            
//            valVsConstIdx_[vals.get(0).intValue()] = idx2D;
//            
//            bPlaceFound = true;
//        }
//
//        if(!bPlaceFound){
//            for (int i = 0; i < chromeConstraints.size(); i++) { //maxslots
//                tempAL = chromeConstraints.get(i);            
//                totalConflictStudents = 0;
//                for (int j = 0; j < tempAL.size(); j++) {
//                    totalConflictStudents += conflictMatrix_[course][tempAL.get(j)];                
//                }              
//                conflictStudents.add(totalConflictStudents);            
//            }           
//            
//            for (int i = 0; i < chromeConstraints.size(); i++) {
//                Prefs.add(Integer.MAX_VALUE);
//            }
//            
//            for (int i : MyRandom.randperm(0,chromeConstraints.size()-1)) {
//            //check surrounding
//            //for (int wt = 0; wt < chromeConstraints.size(); wt++) {
//                worstPref = 0; //best one
//                for (int j = i - MAX_DIST; j <= i+MAX_DIST; j++) {
//                    if(j<0 || j > chromeConstraints.size()-1){
//                        continue;
//                    }
//                    if(conflictStudents.get(j)>0){ //if there is conflict
//                        tempDist = Math.abs(i-j);
//                        if(tempDist == 0){ //&& conflictStudents.get(j) > 0 ... see above
//                            worstPref = Integer.MAX_VALUE;
//                            hardConstraintViolated++;
//                            break;
//                        }
//                        
//                        if(userInput_.bWeighted){ // fitness based on penalty weight function
//                            worstPref+=penaltyFn(tempDist)*conflictStudents.get(j);
////                            if(penaltyFn(tempDist)*conflictStudents.get(j)>worstPref)
////                                worstPref = penaltyFn(tempDist)*conflictStudents.get(j);
//                        }else{ // fitness based on preference                    
//                            if(prefFn(tempDist)>worstPref)
//                                worstPref = prefFn(tempDist);
//                        }   
//                    }
//                }
//                
//                //itry placing Integer.MAX_VALUE in place of getCurPref() below....
//                
//                if(worstPref <= getCurPref()){
//                    chromeConstraints.get(i).add(course);
//                    idx2D = new Idx2D();
//                    idx2D.col = i;
//                    idx2D.position = chromeConstraints.get(i).size()-1;
//                    valVsConstIdx_[course] = idx2D;
//                    bPlaceFound = true;
//                    break;
//                }else{
//                    if(worstPref<this.nextPrefLimit){
//                        this.nextPrefLimit = worstPref;
//                    }
//                }
//            }
//        }
//        
////        int sz = 0;        
////        for (ArrayList consGroup : chromeConstraints) {
////            sz+= consGroup.size();
////        }
//                
//        
//        if(bPlaceFound){
//            
////            if(Prefs.isEmpty())
////                worstPref = 0;
////            else
////                worstPref = Prefs.get(minIdx);
//            curOjbectiveFn(MAX_PREF_WT,worstPref, coeffConflictVsCourse[course], fitness);   
//            //debug(chromeConstraints, fitness);
//            noGood.remove(course*1.0);
//        }else{ 
//            if(hardConstraintViolated == chromeConstraints.size()){//NOGOOD solution eee 
//                noGood.add(course*1.0);//last element
//                
//                if(noGood.subList(0, noGood.size()-1).contains(course*1.0)){
//                    noGood.remove(noGood.size() - 1);//last element
//                    //fitness and satisfactions not changed.
//                }
//            }
//            vals.remove(vals.size() - 1);//last element
//        }        
//    }
//
//   
//    
//    @Override
//    protected boolean isHighlyConstrained(Object obj) {
//        int course = (Integer)obj;
//        
//        if(coeffConflictVsCourse[course]>=userInput_.totalConstraints-maxSlots/4){//top 5
//            return true;
//        }else
//            return false;
//    }
//
//    
//    
//    @Override
//    protected boolean isViolated(final Object obj1, final Object obj2, final Object... additionalInfo) {
//        return hardConstraintViolated((Integer)obj1, (Integer)obj2);
//    }
//    
//    
//    
//    private boolean hardConstraintViolated(final int course1, final int course2){
//        if(conflictMatrix_[course1][course2]==0){
//            return false;
//        }else{
//            return true;
//        }
//    }
//          
//    private void debug(final ArrayList<ArrayList> sat, final ArrayList<Double> fitness){
//        
//////        ArrayList<ArrayList<Integer>> sat = new ArrayList<ArrayList<Integer>>();
//////        for (int wt = 0; wt < 23; wt++) {
//////            sat.add(new ArrayList<Integer>());            
//////        }
//////        
//////        Integer a0[] = {214,	216,	215,	84,	143,	107,	164,	162,	33,	254,	89,	59,	112,	10,	166,	238,	259,	136,	178,	67,	60,	186};
//////        ArrayList<Integer> x0 = new ArrayList<Integer>(Arrays.asList(a0));
//////        sat.get(0).addAll(new ArrayList<Integer>(x0));
//////        
//////        
//////        Integer a1[] = {44,	155,	234,	253,	158,	194,	130,	125,	119};
//////        ArrayList<Integer> x1 = new ArrayList<Integer>(Arrays.asList(a1));
//////        sat.get(1).addAll(new ArrayList<Integer>(x1));
//////        
//////        Integer a2[] = {140,	14,	81,	248,	224,	206,	78,	70,	5,	236};
//////        ArrayList<Integer> x2 = new ArrayList<Integer>(Arrays.asList(a2));
//////        sat.get(2).addAll(new ArrayList<Integer>(x2));  
//////
//////        Integer a3[] = {103,	188,	92,	233,	102,	170,	261};
//////        ArrayList<Integer> x3 = new ArrayList<Integer>(Arrays.asList(a3));
//////        sat.get(3).addAll(new ArrayList<Integer>(x3)); 
//////        
//////        Integer a4[] = {203,	111,	159,	30,	221,	180,	80,	45,	35,	16,	146};
//////        ArrayList<Integer> x4 = new ArrayList<Integer>(Arrays.asList(a4));
//////        sat.get(4).addAll(new ArrayList<Integer>(x4)); 
//////        
//////        Integer a5[] = {83,	46,	28,	129,	126,	22,	252,	207,	240};
//////        ArrayList<Integer> x5 = new ArrayList<Integer>(Arrays.asList(a5));
//////        sat.get(5).addAll(new ArrayList<Integer>(x5)); 
//////        
//////        Integer a6[] = {156,	154,	189,	27,	71,	104,	7,	79,	198,	228,	56,	26};
//////        ArrayList<Integer> x6 = new ArrayList<Integer>(Arrays.asList(a6));
//////        sat.get(6).addAll(new ArrayList<Integer>(x6)); 
//////        
//////        Integer a7[] = {123,	15,	219,	153,	161,	86,	63,	147};
//////        ArrayList<Integer> x7 = new ArrayList<Integer>(Arrays.asList(a7));
//////        sat.get(7).addAll(new ArrayList<Integer>(x7)); 
//////        
//////        Integer a8[] = {256,	110,	230,	176,	157,	47,	205,	171,	213};
//////        ArrayList<Integer> x8 = new ArrayList<Integer>(Arrays.asList(a8));
//////        sat.get(8).addAll(new ArrayList<Integer>(x8)); 
//////        
//////        Integer a9[] = {168,	193,	192,	100,	150,	211,	31,	225,	131,	62};
//////        ArrayList<Integer> x9 = new ArrayList<Integer>(Arrays.asList(a9));
//////        sat.get(9).addAll(new ArrayList<Integer>(x9)); 
//////        
//////        Integer a10[] = {19,	184,	174,	182,	2,	173,	106,	141,	87};
//////        ArrayList<Integer> x10 = new ArrayList<Integer>(Arrays.asList(a10));
//////        sat.get(10).addAll(new ArrayList<Integer>(x10)); 
//////        
//////        Integer a11[] = {220,	50,	113,	43,	54,	172,	200,	241,	137};
//////        ArrayList<Integer> x11 = new ArrayList<Integer>(Arrays.asList(a11));
//////        sat.get(11).addAll(new ArrayList<Integer>(x11)); 
//////        
//////        Integer a12[] = {101,	210,	51,	127,	116,	120,	250,	132,	52,	75,	32,	90,	58,	17};
//////        ArrayList<Integer> x12 = new ArrayList<Integer>(Arrays.asList(a12));
//////        sat.get(12).addAll(new ArrayList<Integer>(x12)); 
//////        
//////        Integer a13[] = {11,	167,	227,	165,	209,	199,	64,	197,	48};
//////        ArrayList<Integer> x13 = new ArrayList<Integer>(Arrays.asList(a13));
//////        sat.get(13).addAll(new ArrayList<Integer>(x13)); 
//////        
//////        Integer a14[] = {12,	144,	82,	29,	212,	1,	222,	4,	122};
//////        ArrayList<Integer> x14 = new ArrayList<Integer>(Arrays.asList(a14));
//////        sat.get(14).addAll(new ArrayList<Integer>(x14)); 
//////        
//////        Integer a15[] = {109,	115,	175,	231,	183,	39,	242};
//////        ArrayList<Integer> x15 = new ArrayList<Integer>(Arrays.asList(a15));
//////        sat.get(15).addAll(new ArrayList<Integer>(x15)); 
//////        
//////        Integer a16[] = {96,	95,	163,	191,	117,	74,	255,	88,	202,	195,	128,	232,	9,	21,	244,	24};
//////        ArrayList<Integer> x16 = new ArrayList<Integer>(Arrays.asList(a16));
//////        sat.get(16).addAll(new ArrayList<Integer>(x16)); 
//////        
//////        Integer a17[] = {94,	6,	114,	105,	133,	217,	258,	142,	229,	249,	135,	118,	148,	34};
//////        ArrayList<Integer> x17 = new ArrayList<Integer>(Arrays.asList(a17));
//////        sat.get(17).addAll(new ArrayList<Integer>(x17)); 
//////        
//////        Integer a18[] = {72,	99,	204,	251,	41};
//////        ArrayList<Integer> x18 = new ArrayList<Integer>(Arrays.asList(a18));
//////        sat.get(18).addAll(new ArrayList<Integer>(x18)); 
//////        
//////        Integer a19[] = {260,	139,	13,	98,	73,	190,	69,	223,	85,	237,	36,	243};
//////        ArrayList<Integer> x19 = new ArrayList<Integer>(Arrays.asList(a19));
//////        sat.get(19).addAll(new ArrayList<Integer>(x19)); 
//////        
//////        Integer a20[] = {245,	196,	151,	160,	187,	55,	49,	257,	97,	65,	61};
//////        ArrayList<Integer> x20 = new ArrayList<Integer>(Arrays.asList(a20));
//////        sat.get(20).addAll(new ArrayList<Integer>(x20)); 
//////        
//////        Integer a21[] = {38,	246,	93,	226,	68,	25,	124,	76,	208,	235,	40,	37};
//////        ArrayList<Integer> x21 = new ArrayList<Integer>(Arrays.asList(a21));
//////        sat.get(21).addAll(new ArrayList<Integer>(x21)); 
//////        
//////        Integer a22[] = {108,	247,	218,	201,	53,	121,	42,	3,	169,	66,	57,	8,	18,	23,	149,	145,	239,	77,	185,	134,	177,	91,	152,	20,	138,	179,	181};
//////        ArrayList<Integer> x22 = new ArrayList<Integer>(Arrays.asList(a22));
//////        sat.get(22).addAll(new ArrayList<Integer>(x22)); 
//////   
//////        System.out.println("\n"+sat);
//////        
//////        for (ArrayList<Integer> list : sat) {
//////            for (int wt = 0; wt < list.size(); wt++) {
//////                list.set(wt, list.get(wt)-1);                
//////            }
//////        }
//        
//        //System.out.println("\n"+sat);      
//        
//        fitness.clear();
//        
//        final int MAX_DIST = 5;
//        ArrayList<Integer> conflictStudents = new ArrayList<Integer>();
//        ArrayList<Integer> tempAL;
//        ArrayList<Integer> Prefs = new ArrayList<Integer>();
//        int totalConflictStudents;
//        double pref;
//
//        int course_col;
//        final int MAX_PREF_WT = 5; //prefFn(0);
//               
//        int course_pos = -1;
//        pref = 0; 
//        course_col = -1;
//        for (ArrayList<Integer> courseList : sat) {
//            course_col++; //first one is 0.
//            course_pos = -1;
//            for (Integer course : courseList) {
//                course_pos++;
//                conflictStudents = new ArrayList<Integer>();
//                for (int i = 0; i < sat.size(); i++) { //maxslots
//                    tempAL = sat.get(i);            
//                    int j;
//                    
//                    if(i==course_col){
//                        j = course_pos+1;
//                    }else{
//                        j = 0;
//                    }
//                    
//                    totalConflictStudents = 0;
//                    for (; j < tempAL.size(); j++) {
//                        if(course!=tempAL.get(j))
//                            totalConflictStudents += conflictMatrix_[course][tempAL.get(j)];                
//                    }                  
//                    conflictStudents.add(totalConflictStudents);            
//                }
//
//                pref = 0;
//                
//                
//                
//                for (int j = course_col; j <= course_col+MAX_DIST; j++) {
//                    if(j<0 || j > sat.size()-1){
//                        continue;
//                    }
//                    
//                    if(conflictStudents.get(j)>0){ //if there is conflict
//                        if(userInput_.bWeighted){                        
//                                pref+=1.0*penaltyFn(j-course_col)*conflictStudents.get(j);                        
//                        }else{
//                            if(prefFn(j-course_col)>(int)pref)
//                                pref = prefFn(j-course_col);
//                                //can use break here, because the cur pref will be maximum.
//                        }
//                    }
//                }
//                          
//                curOjbectiveFn(MAX_PREF_WT,(int)pref, coeffConflictVsCourse[course], fitness);     
//
//            }            
//        }
//        
//        //System.out.append("***your fitness val:" + fitness.get(0)*1.0/STUDENT_SIZE);
//        
//        //Application.getInstance().exit();
//    }
//    
//    
//    
//    /**
//     * NOTE: It only works if the partial solution/preferences is NOT taken into account.
//     * this is bad..... don't use it.
//     * @param vals
//     * @param fitness
//     * @param chromeConstraints
//     * @param noGood
//     * @param valVsConstIdx 
//     */
//    @Override
//    protected void ObjectiveFnReset(final ArrayList<Double> vals, final ArrayList<Double> fitness, 
//    final ArrayList<ArrayList> chromeConstraints, final ArrayList<Double> noGood, 
//    final Idx2D[] valVsConstIdx) {        
//        //the following works if the partial solution/preferences is not taken into account.
//        ArrayList<Double> stepByStepVal = new ArrayList<Double>();
//        ArrayList<Double> stepByStepFitness = new ArrayList<Double>();
//        ArrayList<ArrayList> stepByStepCon = new ArrayList<ArrayList>();
//        ArrayList<Double> stepByStepNoGood = new ArrayList<Double>();
//        Idx2D[] stepByStepvalVsConstIdx = new Idx2D[userInput_.totalConstraints];
//        
//        for (int i = 0; i < stepByStepvalVsConstIdx.length; i++) {
//            stepByStepvalVsConstIdx[i] = new Idx2D();
//        }
//        
//        HashSet<Double> hashSet = new HashSet<Double>(vals);
//        vals.clear();
//        vals.addAll(hashSet);
//
//        
//        for (int i = 0; i < vals.size(); i++) {
//            stepByStepVal.add(vals.get(i));
//            ObjectiveFnAppend(stepByStepVal, stepByStepFitness, stepByStepCon, stepByStepNoGood, stepByStepvalVsConstIdx);
//        }
//        
//        vals.clear();
//        vals.addAll(stepByStepVal);
//        
//        chromeConstraints.clear();
//        chromeConstraints.addAll(stepByStepCon);
//        
//        fitness.clear();;
//        fitness.addAll(stepByStepFitness);
//        
//        noGood.clear();
//        noGood.addAll(stepByStepNoGood);
//        
//        for (int i = 0; i < stepByStepvalVsConstIdx.length; i++) {
//            valVsConstIdx[i] = (Idx2D)stepByStepvalVsConstIdx[i].clone();
//        }
//    }
//
//    @Override
//    protected void ObjectiveFnRemove(final ArrayList<Double> vals, final ArrayList<Double> fitness_, 
//        final ArrayList<ArrayList> constraints, final Idx2D[] valVsConstIdx, final int idx) 
//    throws Exception{
//
////        int course_col;
////        int course_pos;
////        
////        course_col = idx.get(0);
////        course_pos = idx.get(1);        
////        
////        constraints.get(course_col).remove(course_pos);
////        fitness_.clear();
////        
////        debug(constraints, fitness_);
//        
//        int course_col;
//        int course_pos;
//        Idx2D idx2D;
//        Integer val = vals.get(idx).intValue();
//        
//         //debug(constraints, fitness_); 
//        
//        try {
//           
//            course_col = valVsConstIdx[val].col;
//            course_pos = valVsConstIdx[val].position;
//            vals.remove(idx);
//            valVsConstIdx[val] = new Idx2D(); //col = -1; position = -1;
//            constraints.get(course_col).remove(course_pos);
//            
//            for (int i = course_pos; i < constraints.get(course_col).size(); i++) {
//                valVsConstIdx[(Integer)constraints.get(course_col).get(i)].position--;                
//            }                                  
//            
//            fitness_.clear();
//            
//            debug(constraints, fitness_);
//        } catch (Exception e) {
//            e.printStackTrace();
//            throw new Exception("Removal Failure!");
//        }
//
//    }
//
//    
//    @Override
//    protected int getConstraintID(Double val) {
//        return val.intValue();
//    }
//
//
//
//     /**
//     * Clone defined for ExternalData is ONLY SHALLOW CLONE.
//     * @return Object.clone();
//     * @throws CloneNotSupportedException
//     */
//    @Override
//    public Object clone() throws CloneNotSupportedException {
//        return super.clone(); //Not making copy of file data
//    }  
//}
