/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package csp;

import java.util.ArrayList;

/**
 *
 * @author s425770
 */
public class MyAlgorithms{
    public static void virusProliferationAlg(double[] existingVirus, double virus, double[] range){
        
    }
    /**
     * Get directions from one point to another in terms of integers.<br>
     * +1 means towards right. <br>
     * -1 means towards left.<br>
     * 0 means in the same place <br>
     * @param pointFrom - the point from where the direction of ponitTo is to be located
     * @param pointTo - the destination point.
     * @return returns the direction from a set of {-1,0,+1}. <br>
     * +1 means towards right. <br>
     * -1 means towards left. <br>
     * 0 means in the same place <br>
     */
    public static ArrayList<Double> getDirection(final ArrayList<Double> pointFrom, final ArrayList<Double> pointTo){
        ArrayList<Double> direction;
        
        if (pointFrom.size() != pointTo.size()){
            throw new UnsupportedOperationException("Both points should have same dimension");
        }
        
        direction = new ArrayList<Double>(pointFrom.size());
        
        for (int i = 0; i < pointFrom.size(); i++) {
            if(pointTo.get(i) >pointFrom.get(i)){
                direction.add(1.0); //[i] = 1;
            }            
            else if (pointTo.get(i) <pointFrom.get(i)) {
                direction.add(-1.0);
            } 
            else {
                direction.add(0.0); //[i] = 0;
            }                               
        }        
        return direction;
    }
    
    public static ArrayList<ArrayList<Double>> surroundingPointsInfo(ArrayList<ArrayList<Double>> existingVirus, Double[] virus, Double[] range){
        int dim = virus.length;
        int inside = 0;
//        boolean hasNeighbor;
        ArrayList<ArrayList<Double>> surroundingVirus = new ArrayList<ArrayList<Double>>();

        if (virus.length != range.length){
            throw new UnsupportedOperationException("virus and range dimesions do not match");
        }

        for (ArrayList<Double> eVirus : existingVirus) {
            inside = 0;
            for (int i = 0; i < dim; i++) { //dimension of each virus
                if((eVirus.get(i) > (virus[i] - 2*range[i])) && (eVirus.get(i) < (virus[i] + 2*range[i])))
                    inside++;
            }

            if (inside == dim){
                surroundingVirus.add(eVirus);
            }
        }

//        if(surroundingVirus.isEmpty())
//            hasNeighbor = false;
//        else
//            hasNeighbor = true;
        
        return surroundingVirus;
    }

    /**
     * surroundingPointsInfo gives information about the surrounding points of
     * a given point/virus
     * @param existingVirus - all existing points in the calculation space
     * @param virus - points of which surrounding information is to be determined
     * @param range - the range on which surrounding points are checked
     * @return returns ArrayList of surrounding points
     */
    public static ArrayList<Chromosome> surroundingVirusInfo(ArrayList<Chromosome> existingVirus, Chromosome virus, Double[] range){
        int dim = virus.getValsCopy().size();
        int inside = 0;
        ArrayList<Chromosome> surroundingVirus = new ArrayList<Chromosome>();

        if (virus.getValsCopy().size() != range.length){
            throw new UnsupportedOperationException("virus and range dimesions do not match");
        }

        for (Chromosome eVirus : existingVirus) {
            inside = 0;
            for (int i = 0; i < dim; i++) { //dimension of each virus
                if(((Double)eVirus.getVals(i) > (virus.getVals(i) - 2*range[i])) && ((Double)eVirus.getVals(i) < (virus.getVals(i) + 2*range[i])))
                    inside++;
            }

            if (inside == dim){
                surroundingVirus.add(eVirus);
            }
        }
        return surroundingVirus;
    }
}
