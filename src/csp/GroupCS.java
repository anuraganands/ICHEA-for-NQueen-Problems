///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package csp;
//
//import java.util.ArrayList;
//import java.util.Collection;
//
///**
// * Group Constraint Satisfaction
// * @author s425770
// */
//public class GroupCS extends ArrayList{
//    private Integer val;
//    private GroupCS(){
//        ;
//    }
//    public GroupCS(int layers){
//        super();
//        this.layers = layers;
//    }
//    private int layers;
//
//    public int getLayers() {
//        return layers;
//    }
//
//    @Override
//    public boolean add(Object e) {    
//        if(layers > 1){
//            return super.add(e);
//        }
//        else if (layers == 1){
//            val = (Integer)e;
//            return true;
//        }
//        else{
//            return false;
//        }
//    }
//
//    @Override
//    public void add(int index, Object element) {
//        throw new UnsupportedOperationException("Does not support this method! Use add(E) instead.");
//    }
//
//    @Override
//    public boolean addAll(Collection c) {
//        throw new UnsupportedOperationException("Does not support this method! Use add(E) instead.");
//    }
//
//    @Override
//    public boolean addAll(int index, Collection c) {
//        throw new UnsupportedOperationException("Does not support this method! Use add(E) instead.");
//    }
//
//    @Override
//    public Object get(int index) {
//        if(layers > 1){
//            return super.get(index);
//        }
//        else if (layers == 1){
//            return val;
//        }
//        else{
//            throw new NoSuchMethodError("Layers should be >= 1");
//        }
//    }
//
//    @Override
//    public String toString() {
//        if(layers > 1){
//            return super.toString();
//        }
//        else if (layers == 1){
//            return val.toString();
//        }
//        else{
//            throw new NoSuchMethodError("Layers should be >= 1");
//        }        
//    }
// 
//}
