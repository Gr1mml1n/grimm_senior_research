package group.grimm_senior_research;

import java.util.*;
import java.lang.Math;

public class body {
    
    //class variables
    ArrayList<Double> current_position, current_velocity;
    ArrayList<Double> sol4_position, sol4_velocity; //these are the fourth order solutions to pos and vel
    ArrayList<Double> sol5_position, sol5_velocity; //these are the fifth order solutions to pos and vel
    ArrayList<Double> k_1, k_2, k_3, k_4, k_5, k_6; //k terms used to calculate velocity in RK algorithms
    ArrayList<Double> j_1, j_2, j_3, j_4, j_5, j_6; //j terms used to calculate position in RK algorithms
    Double mass;
    
    public static double G = 4*Math.pow(Math.PI,2); //this is 4*pi^2 which is G in AU, Msol, and years.
    //note that this means the simulation will always use units of AU, Msol, and years.
    
    
    //constructor
    public body(double ix, double iy, double ivx, double ivy, double m) {
        this.current_position = new ArrayList<Double>();
        this.current_position.add(ix); this.current_position.add(iy);
        this.current_velocity = new ArrayList<Double>();
        this.current_velocity.add(ivx); this.current_velocity.add(ivy);
        this.mass = m;
        
        this.sol4_position = new ArrayList<Double>();
        this.sol4_position.add(0.0); this.sol4_position.add(0.0);
        this.sol4_velocity = new ArrayList<Double>();
        this.sol4_velocity.add(0.0); this.sol4_velocity.add(0.0);
    }
    
    
    
    
    
    //gravity methods --------------------------------------------------------------------------------------------------
    public ArrayList<Double> Gravity(body b, ArrayList<Double> pos) {
        double ag_mag = -(G * b.mass) / Math.pow(vec_Magnitude(vec_Diff(pos, b.current_position)), 2);
        double ag_ang = vec_Angle(vec_Diff(pos, b.current_position));
        return new_Vec(ag_mag, ag_ang);
    }
    
    //alternate method for when we want to calculate gravity at a different position than current_position.
    public ArrayList<Double> Gravity(body b, ArrayList<Double> b_pos, ArrayList<Double> pos) {
        double ag_mag = -(G * b.mass) / Math.pow(vec_Magnitude(vec_Diff(pos, b_pos)), 2);
        double ag_ang = vec_Angle(vec_Diff(pos, b_pos));
        return new_Vec(ag_mag, ag_ang);
    }
    
    //used to calculate gravity in a 3 body scenario.
    public ArrayList<Double> Gravity2(body b_1, body b_2, ArrayList<Double> pos) {
        ArrayList<Double> ag_1 = Gravity(b_1, pos);
        ArrayList<Double> ag_2 = Gravity(b_2, pos);
        return vec_Add(ag_1, ag_2);
    }
   
    public ArrayList<Double> Gravity2(body b_1, ArrayList<Double> b1_pos, body b_2, ArrayList<Double> b2_pos, ArrayList<Double> pos) {
        ArrayList<Double> ag_1 = Gravity(b_1, b1_pos, pos);
        ArrayList<Double> ag_2 = Gravity(b_2, b2_pos, pos);
        return vec_Add(ag_1, ag_2);
    }
    
    
    
    
    
    //runge kutta methods ----------------------------------------------------------------------------------------------
    
    //Runge-Kutta-Fehlberg algorithm that calculates each body's next position and velocity independently.
    public void RK45(body b_1, body b_2, double h) {
        ArrayList<ArrayList<Double>> k_list = new ArrayList<ArrayList<Double>>();
        ArrayList<ArrayList<Double>> j_list = new ArrayList<ArrayList<Double>>();
        
        k_list.add(this.current_velocity);
        this.j_1 = this.current_velocity;
        j_list.add(this.current_position);
        this.k_1 = Gravity2(b_1, b_2, this.current_position);
        
        k_list.add(vec_Scalmult(this.k_1, 1.0/4.0 * h)); //puts k*a into a list
        this.j_2 = vec_Add2(k_list); //since k is an acceleration, k*a gives us velocity. We add that to the current velocity.
        j_list.add(vec_Scalmult(this.j_1, 1.0/4.0 * h)); //we perform a similar method to get position.
        this.k_2 = Gravity2(b_1, b_2, vec_Add2(j_list));
        
        k_list.set(1, vec_Scalmult(this.k_1, 3.0/32.0 * h));
        k_list.add(vec_Scalmult(this.k_2, 9.0/32.0 * h));
        this.j_3 = vec_Add2(k_list);
        j_list.set(1, vec_Scalmult(this.j_1, 3.0/32.0 * h));
        j_list.add(vec_Scalmult(this.j_2, 9.0/32.0 * h));
        this.k_3 = Gravity2(b_1, b_2, vec_Add2(j_list));
        
        k_list.set(1, vec_Scalmult(this.k_1, 1932.0/2197.0 * h));
        k_list.set(2, vec_Scalmult(this.k_2, -7200.0/2197.0 * h));
        k_list.add(vec_Scalmult(this.k_3, 7296.0/2197.0 * h));
        this.j_4 = vec_Add2(k_list);
        j_list.set(1, vec_Scalmult(this.j_1, 1932.0/2197.0 * h));
        j_list.set(2, vec_Scalmult(this.j_2, -7200.0/2197.0 * h));
        j_list.add(vec_Scalmult(this.j_3, 7296.0/2197.0 * h));
        this.k_4 = Gravity2(b_1, b_2, vec_Add2(j_list));
        
        k_list.set(1, vec_Scalmult(this.k_1, 439.0/216.0 * h));
        k_list.set(2, vec_Scalmult(this.k_2, -8.0 * h));
        k_list.set(3, vec_Scalmult(this.k_3, 3680.0/513.0 * h));
        k_list.add(vec_Scalmult(this.k_4, -845.0/4104.0 * h));
        this.j_5 = vec_Add2(k_list);
        j_list.set(1, vec_Scalmult(this.j_1, 439.0/216.0 * h));
        j_list.set(2, vec_Scalmult(this.j_2, -8.0 * h));
        j_list.set(3, vec_Scalmult(this.j_3, 3680.0/513.0 * h));
        j_list.add(vec_Scalmult(this.j_4, -845.0/4104.0 * h));
        this.k_5 = Gravity2(b_1, b_2, vec_Add2(j_list));
        
        k_list.set(1, vec_Scalmult(this.k_1, -8.0/27.0 * h));
        k_list.set(2, vec_Scalmult(this.k_2, 2.0 * h));
        k_list.set(3, vec_Scalmult(this.k_3, -3544.0/2565.0 * h));
        k_list.set(4, vec_Scalmult(this.k_4, 1859.0/4104.0 * h));
        k_list.add(vec_Scalmult(this.k_5, -11.0/40.0 * h));
        this.j_6 = vec_Add2(k_list);
        j_list.set(1, vec_Scalmult(this.j_1, -8.0/27.0 * h));
        j_list.set(2, vec_Scalmult(this.j_2, 2.0 * h));
        j_list.set(3, vec_Scalmult(this.j_3, -3544.0/2565.0 * h));
        j_list.set(4, vec_Scalmult(this.j_4, 1859.0/4104.0 * h));
        j_list.add(vec_Scalmult(this.j_5, -11.0/40.0 * h));
        this.k_6 = Gravity2(b_1, b_2, vec_Add2(j_list));
        
        this.set_sol4(h); //lastly, we set the solutions after calculating all of these terms.
        this.set_sol5(h);
    }
    
    
//  the following code is used to dynamically calculate each RK45 j,k term based on the other bodies' j,k terms.
    public void RK45_jk1(body b_1, body b_2, double h) {
        this.j_1 = this.current_velocity;
        this.k_1 = Gravity2(b_1, b_2, this.current_position);
    }
    
    public void RK45_jk2(body b_1, body b_2, double h) {
        this.j_2 = vec_Add(this.current_velocity,
                vec_Scalmult(this.k_1, 1.0/4.0 * h));
        
        this.k_2 = Gravity2(b_1, vec_Add(b_1.current_position, //adds up each body's j terms to get their temporary position.
                vec_Scalmult(b_1.j_1, 1.0/4.0 * h)),
                
                b_2, vec_Add(b_2.current_position,
                vec_Scalmult(b_2.j_1, 1.0/4.0 * h)),
                
                vec_Add(this.current_position,
                vec_Scalmult(this.j_1, 1.0/4.0 * h)));
    }
    
    public void RK45_jk3(body b_1, body b_2, double h) {
        this.j_3 = vec_Add(this.current_velocity,
                vec_Add(vec_Scalmult(this.k_1, 3.0/32.0 * h),
                vec_Scalmult(this.k_2, 9.0/32.0 * h)));
        
        this.k_3 = Gravity2(b_1, vec_Add( b_1.current_position,
                vec_Add(vec_Scalmult(b_1.j_1, 3.0/32.0 * h),
                vec_Scalmult(b_1.j_2, 9.0/32.0 * h))),
                
                b_2, vec_Add( b_2.current_position,
                vec_Add(vec_Scalmult(b_2.j_1, 3.0/32.0 * h),
                vec_Scalmult(b_2.j_2, 9.0/32.0 * h))),
                
                vec_Add( this.current_position,
                vec_Add(vec_Scalmult(this.j_1, 3.0/32.0 * h),
                vec_Scalmult(this.j_2, 9.0/32.0 * h))));
    }
    
    public void RK45_jk4(body b_1, body b_2, double h) {
        this.j_4 = vec_Add(this.current_velocity,
                vec_Add(vec_Scalmult(this.k_1, 1932.0/2197.0 * h),
                vec_Add(vec_Scalmult(this.k_2, -7200.0/2197.0 * h),
                vec_Scalmult(this.k_3, 7296.0/2197.0 * h))));
        
        this.k_4 = Gravity2(b_1, vec_Add( b_1.current_position,
                vec_Add(vec_Scalmult(b_1.j_1, 1932.0/2197.0 * h),
                vec_Add(vec_Scalmult(b_1.j_2, -7200.0/2197.0 * h),
                vec_Scalmult(b_1.j_3, 7296.0/2197.0 * h)))),
                
                b_2, vec_Add( b_2.current_position,
                vec_Add(vec_Scalmult(b_2.j_1, 1932.0/2197.0 * h),
                vec_Add(vec_Scalmult(b_2.j_2, -7200.0/2197.0 * h),
                vec_Scalmult(b_2.j_3, 7296.0/2197.0 * h)))),
                
                vec_Add( this.current_position,
                vec_Add(vec_Scalmult(this.j_1, 1932.0/2197.0 * h),
                vec_Add(vec_Scalmult(this.j_2, -7200.0/2197.0 * h),
                vec_Scalmult(this.j_3, 7296.0/2197.0 * h)))));
    }
    
    public void RK45_jk5(body b_1, body b_2, double h) {
        this.j_5 = vec_Add(this.current_velocity,
                vec_Add(vec_Scalmult(this.k_1, 439.0/216.0 * h),
                vec_Add(vec_Scalmult(this.k_2, -8.0 * h),
                vec_Add(vec_Scalmult(this.k_3, 3680.0/513.0 * h),
                vec_Scalmult(this.k_4, -845.0/4104.0 * h)))));
        
        this.k_5 = Gravity2(b_1, vec_Add( b_1.current_position,
                vec_Add(vec_Scalmult(b_1.j_1, 439.0/216.0 * h),
                vec_Add(vec_Scalmult(b_1.j_2, -8.0 * h),
                vec_Add(vec_Scalmult(b_1.j_3, 3680.0/513.0 * h),
                vec_Scalmult(b_1.j_4, -845.0/4104.0 * h))))),
                
                b_2, vec_Add( b_2.current_position,
                vec_Add(vec_Scalmult(b_2.j_1, 439.0/216.0 * h),
                vec_Add(vec_Scalmult(b_2.j_2, -8.0 * h),
                vec_Add(vec_Scalmult(b_2.j_3, 3680.0/513.0 * h),
                vec_Scalmult(b_2.j_4, -845.0/4104.0 * h))))),
                
                vec_Add( this.current_position,
                vec_Add(vec_Scalmult(this.j_1, 439.0/216.0 * h),
                vec_Add(vec_Scalmult(this.j_2, -8.0 * h),
                vec_Add(vec_Scalmult(this.j_3, 3680.0/513.0 * h),
                vec_Scalmult(this.j_4, -845.0/4104.0 * h))))));
    }
    
    public void RK45_jk6(body b_1, body b_2, double h) {
        this.j_6 = vec_Add(this.current_velocity,
                vec_Add(vec_Scalmult(this.k_1, -8.0/27.0 * h),
                vec_Add(vec_Scalmult(this.k_2, 2.0 * h),
                vec_Add(vec_Scalmult(this.k_3, -3544.0/2565.0 * h),
                vec_Add(vec_Scalmult(this.k_4, 1859.0/4104.0 * h),
                vec_Scalmult(this.k_5, -11.0/40.0 * h))))));
        
        this.k_6 = Gravity2(b_1, vec_Add( b_1.current_position,
                vec_Add(vec_Scalmult(b_1.j_1, -8.0/27.0 * h),
                vec_Add(vec_Scalmult(b_1.j_2, 2.0 * h),
                vec_Add(vec_Scalmult(b_1.j_3, -3544.0/2565.0 * h),
                vec_Add(vec_Scalmult(b_1.j_4, 1859.0/4104.0 * h),
                vec_Scalmult(b_1.j_5, -11.0/40.0 * h)))))),
                
                b_2, vec_Add( b_2.current_position,
                vec_Add(vec_Scalmult(b_2.j_1, -8.0/27.0 * h),
                vec_Add(vec_Scalmult(b_2.j_2, 2.0 * h),
                vec_Add(vec_Scalmult(b_2.j_3, -3544.0/2565.0 * h),
                vec_Add(vec_Scalmult(b_2.j_4, 1859.0/4104.0 * h),
                vec_Scalmult(b_2.j_5, -11.0/40.0 * h)))))),
                
                vec_Add( this.current_position,
                vec_Add(vec_Scalmult(this.j_1, -8.0/27.0 * h),
                vec_Add(vec_Scalmult(this.j_2, 2.0 * h),
                vec_Add(vec_Scalmult(this.j_3, -3544.0/2565.0 * h),
                vec_Add(vec_Scalmult(this.j_4, 1859.0/4104.0 * h),
                vec_Scalmult(this.j_5, -11.0/40.0 * h)))))));
    } 
//
    
    //sets fourth order solution
    public void set_sol4(double h) {
       ArrayList<ArrayList<Double>> vel_sum = new ArrayList<ArrayList<Double>>();
       vel_sum.add(vec_Scalmult(this.k_1, 25.0/216.0));
       vel_sum.add(vec_Scalmult(this.k_3, 1408.0/2565.0));
       vel_sum.add(vec_Scalmult(this.k_4, 2197.0/4104.0));
       vel_sum.add(vec_Scalmult(this.k_5, -1.0/5.0));
       this.sol4_velocity = vec_Add(vec_Scalmult(vec_Add2(vel_sum), h), this.current_velocity);
       
       ArrayList<ArrayList<Double>> pos_sum = new ArrayList<ArrayList<Double>>();
       pos_sum.add(vec_Scalmult(this.j_1, 25.0/216.0));
       pos_sum.add(vec_Scalmult(this.j_3, 1408.0/2565.0));
       pos_sum.add(vec_Scalmult(this.j_4, 2197.0/4104.0));
       pos_sum.add(vec_Scalmult(this.j_5, -1.0/5.0));
       this.sol4_position = vec_Add(vec_Scalmult(vec_Add2(pos_sum), h), this.current_position);
    }
    
    //sets fifth order solution
    public void set_sol5(double h) {
       ArrayList<ArrayList<Double>> vel_sum = new ArrayList<ArrayList<Double>>();
       vel_sum.add(vec_Scalmult(this.k_1, 16.0/135.0));
       vel_sum.add(vec_Scalmult(this.k_3, 6656.0/12825.0));
       vel_sum.add(vec_Scalmult(this.k_4, 28561.0/56430.0));
       vel_sum.add(vec_Scalmult(this.k_5, -9.0/50.0));
       vel_sum.add(vec_Scalmult(this.k_6, 2.0/55.0));
       this.sol5_velocity = vec_Add(vec_Scalmult(vec_Add2(vel_sum), h), this.current_velocity);
       
       ArrayList<ArrayList<Double>> pos_sum = new ArrayList<ArrayList<Double>>();
       pos_sum.add(vec_Scalmult(this.j_1, 16.0/135.0));
       pos_sum.add(vec_Scalmult(this.j_3, 6656.0/12825.0));
       pos_sum.add(vec_Scalmult(this.j_4, 28561.0/56430.0));
       pos_sum.add(vec_Scalmult(this.j_5, -9.0/50.0));
       pos_sum.add(vec_Scalmult(this.j_6, 2.0/55.0));
       this.sol5_position = vec_Add(vec_Scalmult(vec_Add2(pos_sum), h), this.current_position);
    }
    
    
    //this is the classic RK4 algorithm with no dynamic time step.
    public void simpleRK4(body b_1, body b_2, double h) {
        this.j_1 = this.current_velocity;
        this.k_1 = Gravity2(b_1, b_2, this.current_position);
        
        this.j_2 = vec_Add(this.current_velocity, vec_Scalmult(this.k_1, 0.5*h));
        this.k_2 = Gravity2(b_1, b_2, vec_Add(this.current_position, vec_Scalmult(this.j_1, 0.5*h)));
        
        this.j_3 = vec_Add(this.current_velocity, vec_Scalmult(this.k_2, 0.5*h));
        this.k_3 = Gravity2(b_1, b_2, vec_Add(this.current_position, vec_Scalmult(this.j_2, 0.5*h)));
        
        this.j_4 = vec_Add(this.current_velocity, vec_Scalmult(this.k_3, h));
        this.k_4 = Gravity2(b_1, b_2, vec_Add(this.current_position, vec_Scalmult(this.j_3, h))); 
        
        this.set_RK4(h);
    }
    

//  the following code is used to dynamically calculate RK4 using each j,k term based on other body's j,k terms.
    public void RK4_jk1(body b_1, body b_2, double h) {
        this.j_1 = this.current_velocity;
        this.k_1 = Gravity2(b_1, b_2, this.current_position);
    }
    
    public void RK4_jk2(body b_1, body b_2, double h) {
        this.j_2 = vec_Add(this.current_velocity, vec_Scalmult(this.k_1, 0.5*h));
        this.k_2 = Gravity2(b_1, vec_Add(b_1.current_position, vec_Scalmult(b_1.j_1, 0.5*h)),
                b_2, vec_Add(b_1.current_position, vec_Scalmult(b_1.j_1, 0.5*h)),
                vec_Add(this.current_position, vec_Scalmult(this.j_1, 0.5*h)));
    }
    
    public void RK4_jk3(body b_1, body b_2, double h) {
        this.j_3 = vec_Add(this.current_velocity, vec_Scalmult(this.k_2, 0.5*h));
        this.k_3 = Gravity2(b_1, vec_Add(b_1.current_position, vec_Scalmult(b_1.j_2, 0.5*h)),
                b_2, vec_Add(b_1.current_position, vec_Scalmult(b_1.j_2, 0.5*h)),
                vec_Add(this.current_position, vec_Scalmult(this.j_2, 0.5*h)));
    }
    
    public void RK4_jk4(body b_1, body b_2, double h) {
        this.j_4 = vec_Add(this.current_velocity, vec_Scalmult(this.k_3, h));
        this.k_4 = Gravity2(b_1, vec_Add(b_1.current_position, vec_Scalmult(b_1.j_3, 0.5*h)),
                b_2, vec_Add(b_1.current_position, vec_Scalmult(b_1.j_3, 0.5*h)),
                vec_Add(this.current_position, vec_Scalmult(this.j_3, h)));
    }
//
    
    //sets fourth order solution using the classic RK4 algorithm.
    public void set_RK4(double h) {
       ArrayList<ArrayList<Double>> vel_sum = new ArrayList<ArrayList<Double>>();
       vel_sum.add(vec_Scalmult(this.k_1, 1.0/6.0));
       vel_sum.add(vec_Scalmult(this.k_2, 2.0/6.0));
       vel_sum.add(vec_Scalmult(this.k_3, 2.0/6.0));
       vel_sum.add(vec_Scalmult(this.k_4, 1.0/6.0));
       this.sol4_velocity = vec_Add(vec_Scalmult(vec_Add2(vel_sum), h), this.current_velocity);
       
       ArrayList<ArrayList<Double>> pos_sum = new ArrayList<ArrayList<Double>>();
       pos_sum.add(vec_Scalmult(this.j_1, 1.0/6.0));
       pos_sum.add(vec_Scalmult(this.j_2, 2.0/6.0));
       pos_sum.add(vec_Scalmult(this.j_3, 2.0/6.0));
       pos_sum.add(vec_Scalmult(this.j_4, 1.0/6.0));
       this.sol4_position = vec_Add(vec_Scalmult(vec_Add2(pos_sum), h), this.current_position);
    }
   
    //method to advance the body forward one time step.
    public void advance(double h) {
       this.current_velocity = this.sol4_velocity;
       this.current_position = this.sol4_position;
    }  
    
    
    
    
    
    //Vector Operations---------------------------------------------------------------------------------------------------
    public static double vec_Magnitude(ArrayList<Double> V) { //must be a vector of size 2.
        return Math.hypot(V.get(0), V.get(1));
    }

    public static double vec_Angle(ArrayList<Double> V) { //must be a vector of size 2.
        return Math.atan2(V.get(1), V.get(0));
    }

    public static ArrayList<Double> new_Vec(double mag, double angle) { //creates a vector given a magnitude and angle.
        ArrayList<Double> V = new ArrayList<Double>();
        V.add(mag * Math.cos(angle)); //sets x direction
        V.add(mag * Math.sin(angle)); //sets y direction
        return V;
    }
    
    public static ArrayList<Double> vec_Diff(ArrayList<Double> V_1, ArrayList<Double> V_2) { //returns the vector V_1 - V_2
        ArrayList<Double> result = new ArrayList<Double>();
        result.add((V_1.get(0) - V_2.get(0)));
        result.add((V_1.get(1) - V_2.get(1)));
        return result;
    }
    
    public static ArrayList<Double> vec_Add(ArrayList<Double> V_1, ArrayList<Double> V_2) { //returns the vector V_1 + V_2
        ArrayList<Double> result = new ArrayList<Double>();
        result.add((V_1.get(0) + V_2.get(0)));
        result.add((V_1.get(1) + V_2.get(1)));
        return result;
    }
    
    public static ArrayList<Double> vec_Add2(ArrayList<ArrayList<Double>> A) { //returns the addition of several vectors
        ArrayList<Double> sum = new ArrayList<Double>();
        sum.add(0,0.0); sum.add(1,0.0);
        for(int i = 0; i < A.size(); i++) {
            sum = vec_Add(sum, A.get(i));
        }
        return sum;
    }
    
    public static ArrayList<Double> vec_Scalmult(ArrayList<Double> V, double scalar) { //multiplies the vector by a scalar
        ArrayList<Double> result = new ArrayList<Double>();
        result.add((scalar * V.get(0)));
        result.add((scalar * V.get(1)));
        return result;
    }
    
    public static double vec_Dot(ArrayList<Double> V_1, ArrayList<Double> V_2) {
        return V_1.get(0)*V_2.get(0) + V_1.get(1)*V_2.get(1);
    }
    
    public static double vec_Crossmag(ArrayList<Double> V_1, ArrayList<Double> V_2) {
        return V_1.get(0)*V_2.get(1) - V_1.get(1)*V_2.get(0);
    }
    
}
