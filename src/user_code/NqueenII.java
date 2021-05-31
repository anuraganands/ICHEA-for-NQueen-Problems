/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package user_code;

import csp.Chromosome;
import csp.CspProcess;
import csp.ExternalData;
//import csp.GroupCS;
import csp.Idx2D;
import csp.MyException;
import csp.MyMath;
import csp.SolutionFoundException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.Set;
import java.util.StringTokenizer;
import javax.swing.JOptionPane;

/**
 *
 * @author s425770
 */
public class NqueenII extends ExternalData{
    private int chessBoardSize_;
    
    public NqueenII(String fileName, int populaion, int generation, int curPref, 
            int prevPref, boolean saveChromes, int solutionBy, Class t) 
            throws InstantiationException, IllegalAccessException, MyException {
        super(fileName, populaion, generation, curPref, prevPref, saveChromes, solutionBy, t);       
        this.chessBoardSize_ = 0;
        this.readData();
    }
       

    @Override
    protected int getConstraintID(Double val) {
        return val.intValue();
    }

    @Override
    protected ArrayList<Chromosome> initializeExternalChrmosomes(final int population) throws SolutionFoundException{
        if(userInput_ == null)
            throw new UnsupportedOperationException("User input not initialized");

        Chromosome chrome;
        ArrayList<Integer> randChromes;
        
        ArrayList<Chromosome> newChromes = new ArrayList<Chromosome>();   
        
        //////chromosome_ = new ArrayList<Chromosome>();
        
        //get chromosomes from file.
        //<<
        try {
            File file = new File(new File(".").getCanonicalPath() + "/chromosomes.txt");
            Scanner chromeFile = new Scanner(file);
                
            int TOTAL_CHROMES = chromeFile.nextInt();
            chromeFile.nextLine(); //ignore comments

            int CHROMES_LENGTH = chromeFile.nextInt();
            chromeFile.nextLine(); //ignore

            int sz;
            String tempStr;
            StringTokenizer str;                   
            
            if (CHROMES_LENGTH > userInput_.totalConstraints){
                throw new IOException("stored data has bigger domain!");
            }
            
            while(chromeFile.hasNext()){    
                tempStr = chromeFile.nextLine();
                str = new StringTokenizer(tempStr," ");
                sz = str.countTokens();
                chrome = new Chromosome(userInput_.solutionBy, this);
                for (int i = 0; i < sz; i++){                
                    tempStr = str.nextElement().toString();  
                                   
                    chrome.appendVal(Double.valueOf(tempStr), Double.MAX_VALUE);
                }
                newChromes.add(chrome);       
            } 
            
        } catch (IOException ioe) {
            //System.out.println("Chromosomes NOT generated from the file");
            newChromes = new ArrayList<Chromosome>();
        }catch (NoSuchElementException nsee){
            //System.out.println("Chromosomes NOT generated from the file. It seems the file is empty.");
            newChromes = new ArrayList<Chromosome>();
        }
        //>>
        
        //System.out.println(chromosome_);
        
        for (int i = 0; i < population; i++) {             
            chrome = new Chromosome(userInput_.solutionBy, this);
            chrome.appendVal(initializeCounter, Double.MAX_VALUE);
            newChromes.add(chrome);
            initializeCounter = (initializeCounter+1)%userInput_.totalConstraints;
        }
  
        //userInput_.population = chromosome_.size();
        return newChromes;
    }

//    @Override
//    public ArrayList<Double> negateVal(ArrayList<Double> vals) {
//        ArrayList<Double> negVals = new ArrayList<Double>();        
//        for (int i = 0; i < vals.size(); i++) {
//           negVals.add(userInput_.maxVals.get(i) - vals.get(i));
//        }
//    }

    @Override
    protected void ObjectiveFnRefresh(ArrayList<ArrayList> chromeConstraints, 
    ArrayList<Double> fitness, ArrayList<Double> vals, Idx2D[] valVsConstIdx_) 
    throws SolutionFoundException{
        ObjectiveFnReset(vals, fitness, chromeConstraints, null, valVsConstIdx_);
    }

    @Override
    protected double fitnessValWeightBased(int penalty, ArrayList<Double> fitness) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean getForcedCSPsol(ArrayList<ArrayList> chromeConstraints, boolean bShowProgress) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected boolean isHighlyConstrained(Object obj) {
        return false;
    }

    @Override
    protected int maxPref() {
        return 1;
    }

    @Override
    public void tryForcedCSPsolUpdate(ArrayList<Double> vals, ArrayList<Double> fitness, ArrayList<ArrayList> chromeConstraints, ArrayList<Double> noGood, Idx2D[] valVsConstIdx_, boolean bShowProgress) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
       
    
    
    
    
    
//    @Override
    protected boolean isViolated(Object obj1, Object obj2, Object... additionalInfo) {
        boolean violated; // 1 = > true; 0 = false;

        if (obj1 instanceof Integer && obj2 instanceof Integer){
            ;
        }else{
            throw new ClassCastException("Expecting Integers");
        }
        
        Integer pos1 = (Integer)obj1;
        Integer pos2 = (Integer)obj2;  
        Integer dist = (Integer)additionalInfo[0];
        
        if(pos1 == pos2){ 
            violated = true; //note..............
            return violated;
        }
        
        if(pos1 == pos2){ //same row
            violated = true; //note..............        
        }else if (Math.abs(pos1-pos2) == dist){//diagonally same
            violated = true;
        }else{
            violated = false;
        }
        
        return violated;
    }

    @Override
    protected void readData() throws MyException {
        String val;
        
        val = JOptionPane.showInputDialog(null, "Enter queen problem size (n)", "Nqueen Problem", JOptionPane.QUESTION_MESSAGE);
        if(val == null)
            throw new UnsupportedOperationException("Cannot accept null value");
        chessBoardSize_ = Integer.parseInt(val);
                    
        //set user input values...        
//        userInput_.fileData = true;//???
        userInput_.totalConstraints = chessBoardSize_;
        userInput_.totalDecisionVars = userInput_.totalConstraints; //currently i take it one for all ordinal values..
        userInput_.totalObjectives = 1; //currently i take it one for all ordinal values..
        userInput_.total__updatedConstraints = userInput_.totalConstraints;
        ExternalData.curHstConstNo = userInput_.total__updatedConstraints;

        for (int i = 0; i < userInput_.totalDecisionVars; i++) {
            userInput_.minVals.add(0.0);
            userInput_.maxVals.add(userInput_.totalConstraints*1.0-1);            
        }               
        
        userInput_.validateData();
        //userInput_.population = Math.min(chessBoardSize_*chessBoardSize_, userInput_.population); //it is minimum.

        userInput_.domainVals = new  ArrayList<ArrayList<Double>>();
        for (int i = 0; i < userInput_.totalDecisionVars; i++) {
            userInput_.domainVals.add(new ArrayList<Double>());
            for (Double j = userInput_.minVals.get(i); j <= userInput_.maxVals.get(i); j++) {
                userInput_.domainVals.get(i).add(j);
            }
        }

        userInput_.doMutation = true;
    }

    @Override
    protected void ObjectiveFnRemove(final ArrayList<Double> vals, final ArrayList<Double> fitness_, 
    final ArrayList<ArrayList> constraints, final Idx2D[] valVsConstIdx, final int idx)
    throws Exception{
        constraints.get(idx).clear();
    }
 
    
    @Override
    protected void ObjectiveFnAppend(final ArrayList<Double> vals, final ArrayList<Double> fitness, 
    final ArrayList<ArrayList> chromeConstraints, final ArrayList<Double> noGood,
    final Idx2D[] valVsConstIdx_) throws SolutionFoundException{
        //??ObjectiveFnReset(vals, fitness, chromeConstraints, noGood, valVsConstIdx_);
        
        ArrayList<Integer> gcs;
        final int queenPos = vals.get(vals.size()-1).intValue();
        
        if(vals.size() == 0)
            return;
        
        //remove duplicates...
        if(vals.subList(0, vals.size()-1).contains(queenPos*1.0)){
            vals.remove(vals.size() - 1);//last element
            //fitness and satisfactions not changed.
            return;
        }
        
        for (int i = 0; i < vals.size()-1; i++) {//all but last one
            if(isViolated(queenPos, vals.get(i).intValue(),Math.abs(vals.size()-1-i))){
                vals.remove(vals.size() - 1);//last element
                return;
            }
        }
        
        while(chromeConstraints.size()<userInput_.totalConstraints){
            chromeConstraints.add(new ArrayList<Integer>());
        } 
        while(fitness.size()<2){
            fitness.add(0.0);//to make it of size 2
        }
        int p = 0;
        while(!chromeConstraints.get(p).isEmpty()){
            p++;
        }
        chromeConstraints.get(p).add(queenPos);


        //<<cater dynamic constraints........................
        double hammingDist = -1;        
        ArrayList<Double> intChromCons = new ArrayList<Double>();
        ArrayList<Double> tempDynamicConst;
        
        for (int i = 0; i < userInput_.totalConstraints; i++) {
            if(!chromeConstraints.get(i).isEmpty())
                intChromCons.add(Double.parseDouble(chromeConstraints.get(i).get(0).toString()));
        }
        
        //<<cater for dynamic constraints....
        int currentSize = chromeConstraints.size();
        for (int i = currentSize; i < userInput_.total__updatedConstraints; i++) {
            chromeConstraints.add(new ArrayList<Integer>());            
        }
        //>>        
        ///Collections.sort(strChromCons);           nooooooo......
        
        int curDynamicLoc = -1;
        //check with dynamic tabu constraints.
        final int maxHamDist = (int)(userInput_.totalConstraints/2);
        //intChromCons = new ArrayList<Double>(intChromCons.subList(fitness.get(1).intValue()-1, intChromCons.size()));
        
        for (int i = 0; i < Math.min(maxDynamicConstraints,CspProcess.dynamicConstraints.size()); i++) {            
            ArrayList<ArrayList> consList = CspProcess.dynamicConstraints.get(i);
            tempDynamicConst = new ArrayList<Double>();
            for (int j = 0; j < userInput_.totalConstraints; j++) {
                if(!consList.get(j).isEmpty())
                    tempDynamicConst.add(Double.parseDouble(consList.get(j).get(0).toString()));
                
            }
            
//            if(intChromCons.size()>tempDynamicConst.size()){
//                continue;
//            }
            hammingDist = MyMath.norm(intChromCons,tempDynamicConst,MyMath.DIST_HAMMING, maxHamDist);            
            curDynamicLoc = userInput_.totalConstraints + i; //- 1 + (dynConsList.get(i) % maxDynamicConstraints+1);
             
            if(hammingDist<=maxHamDist){ //0.2*userInput_.totalConstraints){//why 20%????
                chromeConstraints.get(curDynamicLoc).clear(); //violated                
            }else{
                chromeConstraints.get(curDynamicLoc).clear();
                chromeConstraints.get(curDynamicLoc).add(dynConsList.get(i)); //satisfied
            }
        }
                        
        fitness.clear();
        int vios = 0;
        int length = 0;
        for (ArrayList grp : chromeConstraints.subList(0, userInput_.totalConstraints)) {
            if(grp.isEmpty()){
                break;
            }
            length++;//actually it is satisfaction
        }
        //hence... convert satisfaction to violation.
        vios = userInput_.totalConstraints - length;
        
        
        if(vios ==  0 && !CspProcess.bInTransition){            
            fitness.add(0.0);
            throw new SolutionFoundException("At least one CSP solution found");
        }
        
        for (ArrayList grp : chromeConstraints.subList(userInput_.totalConstraints, userInput_.total__updatedConstraints)) {
            if(grp.isEmpty()){
                vios+=userInput_.totalConstraints;
            }
        }
        
        fitness.add(vios*1.0);//userInput_.total__updatedConstraints - vios*1.0);
        fitness.add(length*1.0); //this satisfaction used for next calculation.
        //>>.......................................
        
        //csp sol
        if(vios ==  0 && !CspProcess.bInTransition){            
            throw new SolutionFoundException("At least one CSP solution found");
        }    
    }

    @Override
    public void addTabuConstraint(ArrayList<ArrayList> localOptConstraints) {                                
        while(CspProcess.dynamicConstraints.size()>=maxDynamicConstraints){
            CspProcess.dynamicConstraints.remove(0);
            ExternalData.dynConsList.remove(0);
            userInput_.total__updatedConstraints--;
        }
        
        int ID = userInput_.totalConstraints;
        
        if(!ExternalData.dynConsList.isEmpty()){
            ID = ExternalData.dynConsList.get(ExternalData.dynConsList.size()-1)+1;            
        }
        
        CspProcess.dynamicConstraints.add(localOptConstraints);
        ExternalData.dynConsList.add(ID); // i think because it starts from 0
        userInput_.total__updatedConstraints++;
    }
    
    
    
    
//    @Override
//    protected void ObjectiveFnReset(ArrayList<Double> vals, ArrayList<Double> fitness_, ArrayList<ArrayList> constraints) {
    @Override
    protected void ObjectiveFnReset(final ArrayList<Double> vals, final ArrayList<Double> fitness, 
    final ArrayList<ArrayList> chromeConstraints, final ArrayList<Double> noGood, 
    final Idx2D[] valVsConstIdx) throws SolutionFoundException{ 
        //GroupCS gcs;
        ArrayList<Integer> gcs;

        //fitness.clear();
        //fitness.add(Double.NaN);
        //chromeConstraints.clear(); //satisfactions
        
        Set<Double> s = new LinkedHashSet<Double>(vals);
        vals.clear();
        vals.addAll(s);
        
        
        for (int i = 0; i < vals.size(); i++) {
            for (int j = vals.size()-1; j >i; j--) {        
//            for (int j = i+1; j < vals.size(); j++) {
                 if(i!=j){
                    if(isViolated(vals.get(i).intValue(), vals.get(j).intValue(),Math.abs(i-j))){
                        vals.remove(j);
                        j--;
                    }else{
                        ;
                    }
                }
            }
        }
        
        chromeConstraints.clear();//satisfactions                                    
        for (Double v : vals) {
            //gcs = new GroupCS(1);
            gcs = new ArrayList<Integer>();
            gcs.add(v.intValue());
            chromeConstraints.add(gcs);
        } 

        
        //<<cater dynamic constraints........................
        double hammingDist = -1;        
        ArrayList<Double> intChromCons = new ArrayList<Double>();
        ArrayList<Double> intTempChromCons;
        
        for (ArrayList grp : chromeConstraints) {
            intChromCons.add(Double.parseDouble(grp.get(0).toString()));
        }
        
        //<<cater for dynamic constraints....
        int currentSize = chromeConstraints.size();
        for (int i = currentSize; i < userInput_.total__updatedConstraints; i++) {
            chromeConstraints.add(new ArrayList<Integer>());            
        }
        //>>        
        ///Collections.sort(strChromCons);           nooooooo......
        
        int curDynamicLoc = -1;
        //check with dynamic tabu constraints.
        final int maxHamDist = 2;
        for (int i = 0; i < Math.min(maxDynamicConstraints,CspProcess.dynamicConstraints.size()); i++) {            
            ArrayList<ArrayList> consList = CspProcess.dynamicConstraints.get(i);
            intTempChromCons = new ArrayList<Double>();
            for (ArrayList grp : consList) {
                if(!grp.isEmpty())
                    intTempChromCons.add(Double.parseDouble(grp.get(0).toString()));
            }
            hammingDist = MyMath.norm(new ArrayList<Double>(intChromCons.subList(0, Math.min(intChromCons.size(),userInput_.totalConstraints))),
            new ArrayList<Double>(intTempChromCons.subList(0, Math.min(intTempChromCons.size(),userInput_.totalConstraints))), 
            MyMath.DIST_HAMMING, maxHamDist);            
            curDynamicLoc = userInput_.totalConstraints + i; //- 1 + (dynConsList.get(i) % maxDynamicConstraints+1);
             
            if(hammingDist<=maxHamDist){ //0.2*userInput_.totalConstraints){//why 20%????
                chromeConstraints.get(curDynamicLoc).clear(); //violated                
            }else{
                chromeConstraints.get(curDynamicLoc).clear();
                chromeConstraints.get(curDynamicLoc).add(dynConsList.get(i)); //satisfied
            }
        }
                        
        fitness.clear();
        int vios = 0;        
        int length = 0;
        for (ArrayList grp : chromeConstraints.subList(0, userInput_.totalConstraints)) {
            if(grp.isEmpty()){
                break;
            }
            length++;//actually it is satisfaction
        }
        //hence... convert satisfaction to violation.
        vios = userInput_.totalConstraints - length;
        
        
        if(vios ==  0 && !CspProcess.bInTransition){            
            fitness.add(0.0);
            throw new SolutionFoundException("At least one CSP solution found");
        }
        
        for (ArrayList grp : chromeConstraints.subList(userInput_.totalConstraints, userInput_.total__updatedConstraints)) {
            if(grp.isEmpty()){
                vios+=userInput_.totalConstraints;
            }
        }
        
        fitness.add(vios*1.0);//userInput_.total__updatedConstraints - vios*1.0);
        fitness.add(length*1.0); //this satisfaction used for next calculation.
        //>>.......................................
        
        //csp sol
        if(vios ==  0 && !CspProcess.bInTransition){            
            throw new SolutionFoundException("At least one CSP solution found");
        }
    }  
}





//
//
///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package user_code;
//
//import csp.Chromosome;
//import csp.ExternalData;
//import csp.MyException;
//import csp.MyRandom;
//import java.io.File;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.HashSet;
//import java.util.NoSuchElementException;
//import java.util.Scanner;
//import java.util.StringTokenizer;
//import javax.swing.JOptionPane;
//
///**
// *
// * @author s425770
// */
//public class NqueenII extends ExternalData{
//    private int chessBoardSize_;
//    
//    public NqueenII(Scanner dataFile, int populaion, int generation, boolean saveChromes, int solutionBy, Class t) throws InstantiationException, IllegalAccessException, MyException {
//        super(dataFile, populaion, generation, saveChromes, solutionBy, t);
//        this.chessBoardSize_ = 0;
//        this.readData();
//    }
//    
//    @Override
//    protected double degreeOfViolation(Object obj1, Object obj2) {                
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    @Override
//    protected int getConstraintID(Double val) {
//        return val.intValue();
//    }
//
//    @Override
//    protected ArrayList<Chromosome> initializeExternalChrmosomes(final int population) {
//        if(userInput_ == null)
//            throw new UnsupportedOperationException("User input not initialized");
//
//        Chromosome chrome;
//        ArrayList<Integer> randChromes;
//        
//        chromosome_ = new ArrayList<Chromosome>();
//        
//        //get chromosomes from file.
//        //<<
//        try {
//            File file = new File(new File(".").getCanonicalPath() + "/chromosomes.txt");
//            Scanner chromeFile = new Scanner(file);
//                
//            int TOTAL_CHROMES = chromeFile.nextInt();
//            chromeFile.nextLine(); //ignore comments
//
//            int CHROMES_LENGTH = chromeFile.nextInt();
//            chromeFile.nextLine(); //ignore
//
//            int sz;
//            String tempStr;
//            StringTokenizer str;                   
//            
//            if (CHROMES_LENGTH > userInput_.totalConstraints){
//                throw new IOException("stored data has bigger domain!");
//            }
//            
//            while(chromeFile.hasNext()){    
//                tempStr = chromeFile.nextLine();
//                str = new StringTokenizer(tempStr," ");
//                sz = str.countTokens();
//                chrome = new Chromosome(userInput_.solutionBy, this);
//                for (int i = 0; i < sz; i++){                
//                    tempStr = str.nextElement().toString();  
//                                   
//                    chrome.appendVal(Double.valueOf(tempStr));
//                }
//                chromosome_.add(chrome);       
//            } 
//            
//        } catch (IOException ioe) {
//            //System.out.println("Chromosomes NOT generated from the file");
//            chromosome_ = new ArrayList<Chromosome>();
//        }catch (NoSuchElementException nsee){
//            //System.out.println("Chromosomes NOT generated from the file. It seems the file is empty.");
//            chromosome_ = new ArrayList<Chromosome>();
//        }
//        //>>
//        
//        randChromes = MyRandom.randperm(0, userInput_.totalConstraints-1);
//        for (int i = 0; chromosome_.size() < population; i++) {            
//            chrome = new Chromosome(userInput_.solutionBy, this);
//            chrome.appendVal(randChromes.get(i%userInput_.totalConstraints).doubleValue());
//            chromosome_.add(chrome);                     
//        }
//
//        return chromosome_;
//    }     
//    
//    @Override
//    protected int isViolated(Object obj1, Object obj2, Object additionalInfo) {
//        int violated; //0 => false; 1 => true
//
//        if (obj1 instanceof Integer && obj2 instanceof Integer){
//            ;
//        }else{
//            throw new ClassCastException("Expecting Integers");
//        }
//        
//        Integer pos1 = (Integer)obj1;
//        Integer pos2 = (Integer)obj2;  
//        Integer dist = (Integer)additionalInfo;
//        
//        if(pos1 == pos2){ 
//            violated = 1; //note..............
//            return violated;
//        }
//        
//        if(pos1 == pos2){ //same row
//            violated = 1; //note..............        
//        }else if (Math.abs(pos1-pos2) == dist){//diagonally same
//            violated = 1;
//        }else{
//            violated = 0;
//        }
//        
//        return violated;
//    }
//
//    @Override
//    protected void readData() throws MyException {
//        String val;
//        
//        val = JOptionPane.showInputDialog(null, "Enter queen problem size (n)", "Nqueen Problem", JOptionPane.QUESTION_MESSAGE);
//        if(val == null)
//            throw new UnsupportedOperationException("Cannot accept null value");
//        chessBoardSize_ = Integer.parseInt(val);
//                    
//        //set user input values...        
//        userInput_.fileData = true;//???
//        userInput_.totalConstraints = chessBoardSize_;
//        userInput_.totalDecisionVars = chessBoardSize_; //currently i take it one for all ordinal values..
//        userInput_.totalObjectives = 1; //currently i take it one for all ordinal values..
//
//        for (int i = 0; i < userInput_.totalDecisionVars; i++) {
//            userInput_.minVals.add(0.0);
//            userInput_.maxVals.add(userInput_.totalConstraints*1.0-1);            
//        }               
//        
//        userInput_.validateData();
//        //userInput_.population = Math.min(chessBoardSize_*chessBoardSize_, userInput_.population); //it is minimum.
//
//        userInput_.domainVals = new  ArrayList<ArrayList<Double>>();
//        for (int i = 0; i < chessBoardSize_; i++) {
//            userInput_.domainVals.add(new ArrayList<Double>());
//            for (int j = 0; j < chessBoardSize_; j++) {
//                userInput_.domainVals.get(i).add(j*1.0);
//            }
//        }
//
//        userInput_.doMutation = true;
//    }
//    
//}
