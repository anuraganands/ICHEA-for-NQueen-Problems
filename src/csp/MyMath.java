/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package csp;

import java.security.InvalidAlgorithmParameterException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import javax.activation.UnsupportedDataTypeException;
import org.jdesktop.application.Application;

/**
 *
 * @author Anurag
 */
public class MyMath {
    public static final int DIST_EUCLEADIAN = 1;    
    public static final int DIST_HAMMING = 2;
    public static final int DIST_DUPLICATE_MISMATCH = 3;
    private static ArrayList<String> permutes;
    
    /**
     * NOTE: It gives unpredictable results if Double value is very large
     * and cannot be converted to long.
     * @param n
     * @param decimalPlaces
     * @return 
     */
    public static double roundN(Double n, int decimalPlaces) {
        if(decimalPlaces < 0){
            //throw new InvalidAlgorithmParameterException("dcimal Place should be positive");
            System.err.println("decimal Place should be positive");
            Application.getInstance().exit();
        }
        
        String fmt = "#.";
        for (int i = 0; i < decimalPlaces; i++) {
            fmt += "#";
        }
        String num = new DecimalFormat(fmt).format(n);
        return Double.valueOf(num);
        
//        if(n > -1.0E14 && n < 1.0E14){
//            long temp;       
//            temp = Math.round(n*Math.pow(10.0,decimalPlaces));
//            return temp/Math.pow(10.0,decimalPlaces);
//        }else{           
//            String snum = n.toString();
//            int dpIdx;
//            int Eidx;
//            int totalDp;
//
//            dpIdx = snum.indexOf("."); //cannot be < 0 in java
//            Eidx = snum.indexOf("E");
//            String numPart = "";
//            String EPart = ""; 
//
//            if(Eidx<0){
//                totalDp = snum.length()-(dpIdx+1);
//                numPart = snum;
//            }else{
//                totalDp = Eidx - (dpIdx+1);
//                EPart = snum.substring(Eidx+1, snum.length());
//                numPart = snum.substring(0, Eidx);
//            }
//
//            String zeros = "";
//            String round = "";
//            if(decimalPlaces > totalDp){
//                for (int i = 0; i < decimalPlaces-totalDp; i++) {
//                    zeros += "0";
//                }
//
//                numPart += zeros;
//            }else if (decimalPlaces == totalDp){
//                ;//no changes
//            }else{//round
//                round = numPart.substring(dpIdx+decimalPlaces+1, dpIdx+decimalPlaces+2);
//                if(Integer.parseInt(round)>4)
//            }
//        }
    }

    public static ArrayList<String> getBinPermutes(int size){
         permutes = new ArrayList<String>();
        
         printBin("", size);
         return permutes;
    }

    public static void printBin(String soFar, int iterations) {
        if(iterations == 0) {
            permutes.add(soFar);
        }
        else {
            printBin(soFar + "0", iterations - 1);
            printBin(soFar + "1", iterations - 1);
        }
    }
    
    public static ArrayList<ArrayList<Double>> getXrandomBinPermutes(final int size,final int limit){
        ArrayList<ArrayList<Double>> combinations = new ArrayList<ArrayList<Double>>();
        
         //printBin("", size);
        double pmOne; //plus minus one
        //int []bits;
        final int combSize = (int)Math.min(limit, Math.pow(2, size));
        
        for (int i = 0; i < combSize; i++) {
            combinations.add(new ArrayList<Double>());
        }
        
        for (int i = 0; i < combSize; i++) {
            //bits = new int[size];
            for (int j = 0; j < size; j++) {
                pmOne = 1.0;
                if(Math.random() < 0.5){
                    pmOne = -1.0;
                } 
                combinations.get(i).add(pmOne);          
            }            
        }
        
        
        
         return combinations;
    }
    
    public static ArrayList<Double> constMultiplicationToVector(double var, final ArrayList<Double> v){
        ArrayList<Double> result;        
        result = new ArrayList<Double>();
                
        for (int i = 0; i < v.size(); i++) {
            result.add(var*v.get(i));
        }
        
        return result;
    }
    
    
    /**
     * Vector multiplication <br>
     * @param entryByEntry true or false <br>
     * if entryByEntry is true then - entry of v1 is multiplied by same index entry of v2
     * if entryByEntry is false - CURRENTLY NOT IMPLEMENTED
     * @param v1 first vector or scalar
     * @param v2 second vector or scalar
     * @return vector multiplication
     * @throws UnsupportedDataTypeException
     */
    public static ArrayList<Double> vectorMultiplication(boolean entryByEntry, final ArrayList<Double> v1, final ArrayList<Double> v2) throws UnsupportedDataTypeException{
        ArrayList<Double> result;
        int resultSize;  
        
        resultSize = v1.size();
        if (v2.size() > v1.size()){
            resultSize = v2.size();
        }
        
        result = new ArrayList<Double>(resultSize);
        
        if(entryByEntry){
                      
               
            if (v1.size() == 1){
                if (v2.size() == 1){
                    result.add(v1.get(0)*v2.get(0));
                }
                else{
                    for (int i = 0; i < v2.size(); i++) {
                        result.add(v1.get(0)*v2.get(i));
                    }
                }
            }
            else{
                if (v2.size() == 1){
                    for (int i = 0; i < v1.size(); i++) {
                        result.add(v1.get(i) * v2.get(0));
                    }
                }
                else{
                    if(v1.size() != v2.size()){
                        throw new UnsupportedOperationException("Dimensions of both point must be same.");
                    }
                    for (int i = 0; i < v1.size(); i++) {
                        result.add(v1.get(i)*v2.get(i));
                    }
                }
            }
            return result;
        }
        return null;
    }
    
    /**
     * Double norm is similar to Matlab's norm. It calculates distance
     * between 2 n dimensional points/vectors.
     * @param p1 - first n dimensional point
     * @param p2 - second n dimensional point
     * @param distancType Following distance types are accepted: <BR>
     * -DIST_EUCLEADIAN = 1: calculates euclidean distance between two vectors.<br>
     * -DIST_HAMMING = 2: calc hamming distance by checking similarity of each element.
     * <I>For eg.</I><BR>
     * <B> 1 2 3 4 <BR>
     * 1 3 7 <BR>
     * ----------- <BR>
     * 0 1 1 1 = 3<BR>
     * </B>
     * -DIST_DUPLICATE_MISMATCH = 3:<<B> Not working properly></B> Currently only checks
     * the number of mismatches in bigger vector. 
     * For eg. <BR>
     * <B> p1 = {1 2 3 4}, p2 = {1 3 7} <BR>
     * mismatch => {2,4} in p1 because p1 is bigger. norm = 2. <BR></B>
     * @return Euclidean distance between the given 2 points is returned.
     */
    public static double norm(final ArrayList<Double> p1, final ArrayList<Double> p2, final int distanceType, Object ... additionalInfo){
        int size;
        Double sum;
        
        if(distanceType == DIST_EUCLEADIAN){
            if(p1.size() != p2.size()){
                throw new UnsupportedOperationException("Dimensions of both point must be same.");
            }

            size = p1.size();

            sum = 0.0;
            for (int i = 0; i < size; i++) {
                sum += Math.pow(p1.get(i) - p2.get(i),2);
            }
            sum = Math.sqrt(sum);
        }else if(distanceType == DIST_DUPLICATE_MISMATCH){
            int duplicates = 0;
            int min;
            ArrayList<Double> unique = new ArrayList<Double>(p1);
            unique.addAll(p2);
            
            min = p1.size();
            if(p2.size()<p1.size())
                min = p2.size();
            
            HashSet<Double> hashSet = new HashSet<Double>(unique);
            unique = new ArrayList<Double>(hashSet);
            sum = (double)(unique.size() - min);
//            for (Double d : p1) {
//                if(!p2.contains(d)){
//                    duplicates++;
//                }
//            } 
//            for (Double d : p2) {
//                if(!p1.contains(d)){
//                    duplicates++;
//                }
//            } 
//            sum = (double)duplicates;
        }else if(distanceType == DIST_HAMMING){
            int maxV = Math.max(p1.size(), p2.size());
            int minV = Math.min(p1.size(), p2.size());
            int dist = 0;
            int maxDistAllowed = (Integer)additionalInfo[0];
            
            for (int i = 0; i < minV; i++) {
                if(Math.abs(roundN(p2.get(i)-p1.get(i),5))>0){
                    dist++;
                    if(dist>=maxDistAllowed){
                        break;
                    }
                }                
            }
            dist += maxV-minV;
            sum = dist*1.0;
        }else{
            sum = Double.NaN;
        }

        return sum;
    }


  
    
    /**
    * mean - calculates mean of the given ArrayList<Double>
    * @param v - ArrayList<Double> of whom mean is to be found
    * @return return mean/average
    */
    public static double mean(ArrayList v){
        if(v.isEmpty()){
            return Double.NaN;
        }        

        double sum = 0.0;
        for (int i = 0; i < v.size(); i++) {
            sum += Double.parseDouble(v.get(i).toString());
        }
        return sum/v.size();
    }

    /**
     * Double sum does the summation of an array from index minIdx to maxIdx
     * @param obj Array to be summed up
     * @param minIdx - minimum index of the array to be included in the summation
     * @param maxIdx - maximum index of the array to be included in the summation
     * @return summation from minIdx to maxIdx
     */
    public static Double sum(Double[] obj, int minIdx, int maxIdx){
        Double sum = 0.0;
        if (maxIdx > obj.length-1)
            throw new ArrayIndexOutOfBoundsException("maxIdx is greater than array length");
        if (minIdx < 0)
            throw new ArrayIndexOutOfBoundsException("minIdx is less than 0");

        for (int i = minIdx; i <= maxIdx; i++) {
            sum += obj[i];
        }

        return sum;
    }

    /**
     * Same as overloaded method sum(Double[] obj,...) at {@link csp.MyMath}
     * @see csp.MyMath
     * @param obj
     * @param minIdx
     * @param maxIdx
     * @return
     */
     public static Integer sum(Integer[] obj, int minIdx, int maxIdx){
        Integer sum = 0;
        if (maxIdx > obj.length-1)
            throw new ArrayIndexOutOfBoundsException("maxIdx is greater than array length");
        if (minIdx < 0)
            throw new ArrayIndexOutOfBoundsException("minIdx is less than 0");

        for (int i = minIdx; i <= maxIdx; i++) {
            sum += obj[i];
        }

        return sum;
    }   
     
    /**
      * vectorAddition add 2 vectors (one dimesional matrix) of same size
      * @param m1 first vector
      * @param m2 second vector
      * @return summation of 2 vectors
      */ 
    public static ArrayList<Double> vectorAddition(ArrayList<Double> m1, ArrayList<Double> m2){
        ArrayList<Double> result;

        if(m1.size() != m2.size()){
            throw new UnsupportedOperationException("Dimensions of both vectors must be same.");      
        }
        result = new ArrayList<Double>(m1.size());

        for (int i = 0; i < m1.size(); i++) {
            result.add(m1.get(i)+m2.get(i));            
        }
        return result;
    } 
    
    /**
     * Subtracts two vectors element by element (m1 - m2)
     * @param m1 - first vector
     * @param m2 - second vector
     * @return result of (m1-m2)
     */
    public static ArrayList<Double> vectorSubtraction(ArrayList<Double> m1, ArrayList<Double> m2){
        ArrayList<Double> result;

        if(m1.size() != m2.size()){
            throw new UnsupportedOperationException("Dimensions of both vectors must be same.");      
        }
        result = new ArrayList<Double>(m1.size());

        for (int i = 0; i < m1.size(); i++) {
            result.add(m1.get(i)-m2.get(i));            
        }
        return result;
    } 
    
    /**
     * Find <b>Index</b> of min value of any type of given array
     * @param in ArrayList must be instance of Comparable
     * @return return index of first found minimum element
     */
    public static int minIdx(ArrayList<? extends Comparable> in){        
        Comparable minVal;
        int minIdx = 0;
        
        if(in.isEmpty()){
            minVal = -1;
        }        
        
        minVal = in.get(minIdx);
        
        for (int i = 1; i < in.size(); i++) {
            if(minVal.compareTo(in.get(i))>0){
                minVal = in.get(i);
                minIdx = i;
            }            
        }  
        
        return minIdx;
    }
}
