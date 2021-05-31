///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package user_code;
//
//import csp.ExternalData;
//import csp.MyException;
//import java.util.ArrayList;
//import java.util.Scanner;
//import java.util.StringTokenizer;
//
///**
// *
// * @author s425770
// */
//public class PerfectSquare extends ExternalData{
//    private int problemNo_;
//    private int totalSq_;
//    private int masterSq_;
//    ArrayList<Integer> squares_;
//    
//    public PerfectSquare(Scanner dataFile, int populaion, int generation, boolean saveChromes, int solutionBy, Class t) throws InstantiationException, IllegalAccessException, MyException {
//        super(dataFile, populaion, generation, saveChromes, solutionBy, t);
//        this.totalSq_ = -1;
//        this.masterSq_ = -1;
//        squares_ = new ArrayList<Integer>();
//        this.readData();
//    }
//    
//    @Override
//    protected void readData() throws MyException {
//       while (dataFile_.hasNext("#")){                
//            dataFile_.nextLine(); //ignore comments  
//        }
//
//        StringTokenizer str;
//        Integer neighborState;
//        int nextState = 0;
//        int sz;
//        String tempStr;
//    
//        tempStr = dataFile_.nextLine();
//        str = new StringTokenizer(tempStr," ");
//        sz = str.countTokens();
//        if(sz<4){
//            throw new MyException("Incorrect File Data. At least 4 integers are required.");
//        }
//        
//        tempStr = str.nextElement().toString();
//        this.problemNo_ = Integer.valueOf(tempStr);
//        
//        tempStr = str.nextElement().toString();
//        this.totalSq_ = Integer.valueOf(tempStr);
//        
//        tempStr = str.nextElement().toString();
//        this.masterSq_ = Integer.valueOf(tempStr);        
//        
//        for (int i = 0; i < sz-3; i++){
//            tempStr = str.nextElement().toString();
//            squares_.add(Integer.valueOf(tempStr));
//        }       
//
//        userInput_.fileData = true;
//        userInput_.totalConstraints = this.totalSq_;
//        userInput_.totalDecisionVars = userInput_.totalConstraints;
//        userInput_.totalObjectives = 0;
//        for (int i = 0; i < userInput_.totalDecisionVars; i++) {
//            userInput_.minVals.add(1.0);
//            userInput_.maxVals.add(userInput_.totalConstraints*1.0);            
//        }        
//        userInput_.validateData();
//
//        userInput_.domainVals = new  ArrayList<ArrayList<Double>>();
//        for (int i = 1; i <= userInput_.totalConstraints; i++) {
//            userInput_.domainVals.add(new ArrayList<Double>());
//            for (int j = 1; j <= userInput_.totalDecisionVars; j++) {
//                userInput_.domainVals.get(i).add(j*1.0);
//            }
//        }
//
//        userInput_.doMutation = true;
//    }
//    
//    @Override
//    protected int getConstraintID(Double val) {
//        return val.intValue();
//    }
//    
//    @Override
//    protected double degreeOfViolation(Object obj1, Object obj2) {                
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//    
//    @Override
//    protected boolean isViolated(Object obj1, Object obj2, Object additionalInfo) {
//        boolean violated;
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
//            violated = true; //note..............
//            return violated;
//        }
//        
//        if(pos1 == pos2){ //same row
//            violated = true; //note..............        
//        }else if (Math.abs(pos1-pos2) == dist){//diagonally same
//            violated = true;
//        }else{
//            violated = false;
//        }
//        
//        return violated;
//    }
//}
