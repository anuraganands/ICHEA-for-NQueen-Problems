/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package csp;

/**
 *
 * @author s425770
 */
public class SolutionFoundException extends Exception{
    private String message_;
    public SolutionFoundException(String message){
        super();
        message_ = message;
    }
    
    @Override
    public String getMessage() {
        return this.message_;
    }    
    
    public void printMessage(){
        System.out.flush(); // It is must otherwise there is dealy in printing System.err
        System.err.flush();
        System.err.println(this.message_);
        System.err.flush();
        System.out.flush();
    }

    @Override
    public String toString() {
        return this.message_;
    }
}
