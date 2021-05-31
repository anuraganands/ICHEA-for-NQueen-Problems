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
//public class Nqueen extends ExternalData{
//    private int chessBoardSize_;
//    private Integer encodeCapacity_;
//    
//    public Nqueen(Scanner dataFile, int populaion, int generation, boolean saveChromes, int solutionBy, Class t) throws InstantiationException, IllegalAccessException, MyException {
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
//        return decode(val.intValue(), 'r');
//    }
//
//    @Override
//    protected ArrayList<Chromosome> initializeExternalChrmosomes(final int population) {
//        if(userInput_ == null)
//            throw new UnsupportedOperationException("User input not initialized");
//
//        Chromosome chrome;
//        
//        double [][]temp = new double[userInput_.totalConstraints][userInput_.totalConstraints]; 
//        int tempRow, tempCol;
//                
//        for (Integer i = 0; i < userInput_.totalConstraints; i++) {            
//            for (int j = 0; j < userInput_.domainVals.get(i).size(); j++) {
////                chrome = new Chromosome(Chromosome.BY_SATISFACTIONS, this);
////                chrome.appendVal(userInput_.domainVals.get(i).get(j));
////                chromosome_.add(chrome);
//                temp[i][j] = userInput_.domainVals.get(i).get(j);
//            }            
//        }
//        ArrayList<Double> vals = new ArrayList<Double>();
//        int beforeSize, afterSize;
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
//            int length;
//            int r, c;            
//            int oldEncodeCapacity = encodeCapacity_;                         
//            length = Integer.toString(CHROMES_LENGTH).concat(Integer.toString(CHROMES_LENGTH)).length();
//            int newEncodeCapacity = (int)Math.pow(10.0, length);
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
//                    encodeCapacity_ = newEncodeCapacity; //(int)Math.pow(10.0, length);
//                    r = decode(Double.valueOf(tempStr).intValue(), 'r');
//                    c = decode(Double.valueOf(tempStr).intValue(), 'c');
//                    
//                    encodeCapacity_ = oldEncodeCapacity;                    
//                    chrome.appendVal(encode(r, c).doubleValue());
//                }
//                chromosome_.add(chrome);       
//            } 
//            encodeCapacity_ = oldEncodeCapacity;
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
//        //System.out.println(chromosome_);
//        
//        for (int i = 0; chromosome_.size() < population; i++) {
//            tempRow = MyRandom.randperm(0, userInput_.totalConstraints-1).get(0);
//            tempCol = MyRandom.randperm(0, userInput_.totalConstraints-1).get(0);
//            
//            chrome = new Chromosome(userInput_.solutionBy, this);
//            chrome.appendVal(userInput_.domainVals.get(tempRow).get(tempCol));
//                        
//            beforeSize = vals.size();            
//            //vals.add(userInput_.domainVals.get(tempRow).get(tempCol));
//            vals.add((double)tempRow);
//
//            HashSet<Double> hashSet = new HashSet<Double>(vals);
//            vals = new ArrayList<Double>(hashSet);            
//            afterSize = vals.size();
//            
//            if(afterSize>beforeSize || vals.size() >= userInput_.totalConstraints ){
//                chromosome_.add(chrome);
//            }            
//        }
//
//        
//        //userInput_.population = chromosome_.size();
//        return chromosome_;
//    }
//
////    @Override
////    public ArrayList<Double> negateVal(ArrayList<Double> vals) {
////        ArrayList<Double> negVals = new ArrayList<Double>();        
////        for (int i = 0; i < vals.size(); i++) {
////           negVals.add(userInput_.maxVals.get(i) - vals.get(i));
////        }
////    }
//    
//    
//
//    private Integer encode(Integer row, Integer col){
//        int len = encodeCapacity_.toString().length()-1;//even number
//        len = len/2;
//        String str = "1";
//
//        for (int i= 0; i < len - row.toString().length(); i++) {
//            str += "0";
//        }
//        str += row.toString();
//
//        for (int i = 0; i < len - col.toString().length(); i++) {
//            str += "0";
//        }
//        str+= col.toString();
//        
//        return Integer.parseInt(str);
//    }
//
//    private Integer decode(Integer encodedVal, char rowCol){
//        int length;
//        String str = encodedVal.toString();
//        
//        length = encodeCapacity_.toString().length()-1; //should always be even
//        length = length/2;
//        if(rowCol == 'r'){
//            str = str.substring(1, length+1);
//        }else if(rowCol == 'c'){
//            str = str.substring(length+1, encodeCapacity_.toString().length());
//        }else{
//            str = "";
//            throw new IllegalArgumentException("Incorrect value of rowCol parameter");
//        }
//
//        return Integer.parseInt(str);
//    }
//    
//    @Override
//    protected int isViolated(Object obj1, Object obj2, Object... additionalInfo) {
//        int violated;
//
//        if (obj1 instanceof Integer && obj2 instanceof Integer){
//            ;
//        }else{
//            throw new ClassCastException("Expecting Integers");
//        }
//        
//        Integer pos1 = (Integer)obj1;
//        Integer pos2 = (Integer)obj2;
//
//        int r1, c1;
//        int r2, c2;
//        
//        r1 = decode(pos1, 'r');
//        c1 = decode(pos1, 'c');
//        
//        r2 = decode(pos2, 'r');
//        c2 = decode(pos2, 'c');
//        
//        if(pos1 == pos2){ 
//            violated = 1; //note..............
//            return violated;
//        }
//        
//        if(r1 == r2){ //same row
//            violated = 1; //note..............
//        }else if (c1 == c2) {//same column
//            violated = 1;
//        }else if (Math.abs(r1 - r2) == Math.abs(c1 - c2)){//diagonally same
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
//        int length;
//        
//        val = JOptionPane.showInputDialog(null, "Enter queen problem size (n)", "Nqueen Problem", JOptionPane.QUESTION_MESSAGE);
//        if(val == null)
//            throw new UnsupportedOperationException("Cannot accept null value");
//        chessBoardSize_ = Integer.parseInt(val);
//               
//        length = Integer.toString(chessBoardSize_-1).concat(Integer.toString(chessBoardSize_-1)).length();
//        encodeCapacity_ = (int)Math.pow(10.0, length); //make odd so that 1 is always in the front
//        
//        
//        //set user input values...        
//        userInput_.fileData = true;//???
//        userInput_.totalConstraints = chessBoardSize_;
//        userInput_.totalDecisionVars = chessBoardSize_; //currently i take it one for all ordinal values..
//        userInput_.totalObjectives = 1; //currently i take it one for all ordinal values..
//        String str = "1";
//        for (int i = 0; i < 2; i++) {
//            str += Integer.toString(chessBoardSize_-1);
//        }
//
//        for (int i = 0; i < userInput_.totalDecisionVars; i++) {
//            userInput_.minVals.add(encodeCapacity_.doubleValue());
//            userInput_.maxVals.add(Double.parseDouble(str));            
//        }               
//        
//        userInput_.validateData();
//        userInput_.population = Math.min(chessBoardSize_*chessBoardSize_, userInput_.population); //it is minimum.
//
//        userInput_.domainVals = new  ArrayList<ArrayList<Double>>();
//        for (int i = 0; i < chessBoardSize_; i++) {
//            userInput_.domainVals.add(new ArrayList<Double>());
//            for (int j = 0; j < chessBoardSize_; j++) {
//                userInput_.domainVals.get(i).add(encode(i, j).doubleValue());
//            }
//        }
//
//        userInput_.doMutation = true;
//    }
//    
//}
