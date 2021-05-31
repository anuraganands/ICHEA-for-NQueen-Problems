///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package user_code;
//
//import csp.Chromosome;
//import csp.ExternalData;
//import csp.MyException;
//import csp.MyMath;
//import csp.MyRandom;
//import java.util.ArrayList;
//import java.util.LinkedHashSet;
//import java.util.Scanner;
//import java.util.Set;
//import java.util.StringTokenizer;
//
///**
// *
// * @author s425770
// */
//public class GraphColoring extends ExternalData{
//    private ArrayList<ArrayList<Integer>> nodes_;
//    private Double divider_;
//    
//    public GraphColoring(Scanner dataFile, int populaion, int generation, boolean saveChromes, int solutionBy, Class t) throws InstantiationException, IllegalAccessException, MyException {
//        super(dataFile, populaion, generation, saveChromes, solutionBy, t);
//        this.nodes_ = new ArrayList<ArrayList<Integer>>();
//        this.readData();
//    }
//
//    @Override
//    protected void readData() throws MyException {
//        final int TOTAL_NODES;
//        final Integer TOTAL_COLORS;
//        ArrayList<Integer> arrayList;
//        
//        while (dataFile_.hasNext("c")){                
//            dataFile_.nextLine(); //ignore comments  
//        }
//
//        if(dataFile_.hasNext("p")){
//            dataFile_.next();
//            if(dataFile_.hasNext("edge")){
//                dataFile_.next();
//                TOTAL_NODES = dataFile_.nextInt();
//                TOTAL_COLORS = dataFile_.nextInt();
//            }else{
//                throw new MyException("Incorrect File Format!");
//            }
//        }else{
//            throw new MyException("Incorrect File Format!");
//        }
//        
//        for (int i = 0; i < TOTAL_NODES; i++) {
//            arrayList = new ArrayList<Integer>();
//            nodes_.add(arrayList);  
//        }
//        
//        StringTokenizer str;
//        Integer from;
//        Integer to;
//        int sz;
//        String tempStr;
//        
//        dataFile_.nextLine();
//        while(dataFile_.hasNext("e")){    
//            tempStr = dataFile_.nextLine();
//            str = new StringTokenizer(tempStr," ");
//            sz = str.countTokens();
//            
//            if(sz != 3){
//                throw new MyException("Incorrect File Format!");
//            }
//            
//            str.nextElement(); //ignore e
//            from = Integer.parseInt(str.nextElement().toString())-1;
//            to = Integer.parseInt(str.nextElement().toString())-1;
//            nodes_.get(from).add(to); 
//            nodes_.get(to).add(from);
//        }
//
//        userInput_.fileData = true;
//        userInput_.totalConstraints = TOTAL_NODES;
//        userInput_.totalDecisionVars = userInput_.totalConstraints;;
//        userInput_.totalObjectives = 0
//                ;
//        for (int i = 0; i < userInput_.totalDecisionVars; i++) {
//            userInput_.minVals.add(0.0);
//            userInput_.maxVals.add(TOTAL_NODES-1.0);            
//        }        
//        userInput_.validateData();
//        //userInput_.population = TOTAL_STATES;
//
//        userInput_.domainVals = new  ArrayList<ArrayList<Double>>();
//        Integer lastIdx = TOTAL_COLORS - 1;
//        divider_ = Math.pow(10, lastIdx.toString().length());
//        for (int i = 0; i < userInput_.totalConstraints; i++) {
//            userInput_.domainVals.add(new ArrayList<Double>());
//            for (int j = 0; j < TOTAL_COLORS; j++) {
//                userInput_.domainVals.get(i).add(i + j/divider_);
//            }
//        }
//
//        userInput_.doMutation = true;
//    }
//
////    @Override
////    public ArrayList<Double> negateVal(ArrayList<Double> vals) {
////        throw new UnsupportedOperationException("Not supported yet.");
////    }
//
//    
//    
//    
//    @Override
//    protected int isViolated(Object obj1, Object obj2, Object... additionalInfo) {
//        int violated; // 1 = true, 0 = false;
//        int decimalPlace;
//
//        if (obj1 instanceof Double && obj2 instanceof Double){
//            ;
//        }else{
//            throw new ClassCastException("Expecting Double");
//        }                
//        
//        Double nodeColor1 = (Double)obj1;
//        Double nodeColor2 = (Double)obj2;
//        
//        if(nodeColor1.intValue() == nodeColor2.intValue()){
//            violated = 1; //note..............
//            return violated;
//        }
//
//        if(nodes_.get(nodeColor1.intValue()).isEmpty() || nodes_.get(nodeColor2.intValue()).isEmpty()){
//            violated = 0;
//            return violated;
//        }
//        
//        Integer intVal = divider_.intValue();
//        decimalPlace = intVal.toString().length()-1;
//        
//        if(MyMath.roundN(nodeColor1 - nodeColor1.intValue(), decimalPlace) == 
//                MyMath.roundN(nodeColor2 - nodeColor2.intValue(), decimalPlace)){            
//            if(nodes_.get(nodeColor1.intValue()).contains(nodeColor2.intValue()))
//                violated = 1;
//            else
//                violated = 0;
//        }else{
//            violated = 0;    
//        }
//        
//        return violated;
//    }
//
//    @Override
//    protected double degreeOfViolation(Object obj1, Object obj2) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    @Override
//    protected ArrayList<Chromosome> initializeExternalChrmosomes(final int population) {
//        if(userInput_ == null)
//            throw new UnsupportedOperationException("User input not initialized");
//
//        Chromosome chrome;
//                
//        int j;
//        while(chromosome_.size()<userInput_.population){
//            for (Integer i = 0; i < userInput_.totalConstraints; i++) {            
//                    chrome = new Chromosome(userInput_.solutionBy, this);
//                    j = MyRandom.randperm(0, userInput_.domainVals.get(i).size()-1).get(0);
//                    chrome.appendVal(userInput_.domainVals.get(i).get(j));
//                    chromosome_.add(chrome);         
//                    if(chromosome_.size() >= userInput_.population){
//                        break;
//                    }
//            }
//        }
//
//        //userInput_.population = chromosome_.size();
//        return chromosome_;
//    }
//
//    @Override
//    protected int getConstraintID(Double val) {
//        return val.intValue();
//    }
//
//
//
////    @Override
////    protected Set<Double> domainValues(int dimension) {
////        Set<Double> s = new LinkedHashSet<Double>();
////        for (int j = 0; j < userInput_.totalDecisionVars; j++) {
////            s.add(j/10.0);
////        }
////        return s;
////    }
//}
