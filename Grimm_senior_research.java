package group.grimm_senior_research;

import java.io.*;
import java.util.*;
import java.lang.Math;

public class Grimm_senior_research {

    public static void main(String[] args) {
        
        //creates the master csv file to store summary info of each run.
        try {
            //C:\Users\grimm\OneDrive\Documents\COLLEGE CLASSES\Senior Research\senior_research_code\
            FileWriter master = new FileWriter("master.csv");
            master.write("mass,periapsis,eccentricity,semi-major,period,instability,error\n");
            master.close();
        }
        catch (IOException e) {
            System.out.println("an error occurred.");
        }
        
        
        //this is where to set periapsis and eccentricity ranges. Change these to loops to fit your needs.
        double p_array[] = {4.0};
        double e_array[] = {0.4};
        
        //the simulation will be run for each combination of periapsis and eccentricity. ----------------------------------
        for(double periapsis : p_array) {
        for(double eccentricity : e_array) {
        
            
        //initiation step
        //setting time stuff
        double time = 0.0;
        double tolerance = 0.1/(1.496*Math.pow(10, 11)); //the numerator is the tolerance in meters
        double step = 0.01; //this is the initial step. It will be adapted by the algorithm.
        
        //setting star 1 and 2
        body S_1 = new body(1.0,0.0,0.0, Math.PI,1.0); //star 1
        body S_2 = new body(-1.0,0.0,0.0,-Math.PI,1.0); //star 2
        
        //setting planet
        double a = periapsis*(1 + eccentricity)/(1 - Math.pow(eccentricity, 2));
        double period = Math.sqrt(Math.pow(a, 3)/2);
        double p_velocity = Math.sqrt(4*Math.pow(Math.PI,2)*2/a*(1 + eccentricity)/(1 - eccentricity));
        body Planet = new body(periapsis,0.0,0.0,p_velocity,0.001); //planet
        
        //setting arrays
        ArrayList<ArrayList<Double>> Record = new ArrayList<ArrayList<Double>>(); 
        //Record is a 2d array containing the information at each time step. It is outputted to each run's csv file.
        ArrayList<Double> Row = new ArrayList<Double>();
        //Row is just one row in record.
        
        //fixing and calculating a few values
        com_Vel_Correction(S_1, S_2, Planet);
        double initial_energy = total_Energy(Planet, S_1, S_2);
        double energy_err = 0;
        
        //record initial stuff [time, position, velocity, total energy, planet's momentum, semi-major axis orientation?]
        Row.add(time);
        Row.add(Planet.current_position.get(0)); Row.add(Planet.current_position.get(1));
        Row.add(Planet.current_velocity.get(0)); Row.add(Planet.current_velocity.get(1));
        Row.add(S_1.current_position.get(0)); Row.add(S_1.current_position.get(1));
        Row.add(S_1.current_velocity.get(0)); Row.add(S_1.current_velocity.get(1));
        Row.add(S_2.current_position.get(0)); Row.add(S_2.current_position.get(1));
        Row.add(S_2.current_velocity.get(0)); Row.add(S_2.current_velocity.get(1));
        Row.add(total_Energy(Planet, S_1, S_2));
        Row.add(energy_err);
        Row.add(vec_Crossmag(vec_Diff(Planet.current_position, center_Of_Mass(Planet, S_1, S_2)), Planet.current_velocity)*Planet.mass);
        Record.add((ArrayList<Double>)Row.clone());
        
        
        
        
        //MAIN LOOP TO RUN SIMULATION. ----------------------------------------------------------------------------------
        while(time < 100*period) {
            
            dynamic_RK45(S_1, S_2, Planet, step, tolerance);
            step = Math.min(adapt_Step(S_1, step, tolerance), Math.min(adapt_Step(S_2, step, tolerance), adapt_Step(Planet, step, tolerance)));
            dynamic_RK45(S_1, S_2, Planet, step, tolerance);
            
            
            //uncomment this block to force star 1 and 2 into circular orbits.
//            Planet.RK45(S_1, S_2, step);
//            step = adapt_Step(Planet,step, tolerance);
//         
//            Planet.RK45(S_1, S_2, step); //now that we have the new time step, we find all of the positions.
//            S_1.sol4_position.set(0, Math.cos(Math.PI*(time + step))); S_1.sol4_position.set(1, Math.sin(Math.PI*(time + step)));
//            S_1.sol4_velocity.set(0, -Math.PI*Math.sin(Math.PI*(time + step))); S_1.sol4_velocity.set(1, Math.PI*Math.cos(Math.PI*(time + step)));
//            S_2.sol4_position.set(0, -Math.cos(Math.PI*(time + step))); S_2.sol4_position.set(1, -Math.sin(Math.PI*(time + step)));
//            S_2.sol4_velocity.set(0, Math.PI*Math.sin(Math.PI*(time + step))); S_2.sol4_velocity.set(1, -Math.PI*Math.cos(Math.PI*(time + step)));
            
            //advance forward one step
            Planet.advance(step);
            S_1.advance(step);
            S_2.advance(step);
            com_Vel_Correction(S_1, S_2, Planet);
            
            
            //energy offset for planet (if used)
            energy_err = total_Energy(Planet, S_1, S_2) - initial_energy;
            
            //advance time
            time += step;
            
            //record stuff
            Row.set(0,time);
            Row.set(1,Planet.current_position.get(0)); Row.set(2,Planet.current_position.get(1));
            Row.set(3,Planet.current_velocity.get(0)); Row.set(4,Planet.current_velocity.get(1));
            Row.set(5,S_1.current_position.get(0)); Row.set(6,S_1.current_position.get(1));
            Row.set(7,S_1.current_velocity.get(0)); Row.set(8,S_1.current_velocity.get(1));
            Row.set(9,S_2.current_position.get(0)); Row.set(10,S_2.current_position.get(1));
            Row.set(11,S_2.current_velocity.get(0)); Row.set(12,S_2.current_velocity.get(1));
            Row.set(13,total_Energy(Planet, S_1, S_2));
            Row.set(14, energy_err);
            Row.set(15,vec_Crossmag(vec_Diff(Planet.current_position, center_Of_Mass(Planet, S_1, S_2)), Planet.current_velocity)*Planet.mass);
            Record.add((ArrayList<Double>)Row.clone());
        }
        
        
        
        //write out to run's csv file
        String name = "p-" + String.valueOf(periapsis) + "_e-" + String.valueOf(eccentricity) + ".csv";
        try {
            //C:\Users\grimm\OneDrive\Documents\COLLEGE CLASSES\Senior Research\senior_research_code\
            FileWriter test = new FileWriter(name);
            test.write("time,px,py,pvx,pvy,s1x,s1y,s1vx,s1vy,s2x,s2y,s2vx,s2vy,energy,E_offset,momentum\n"); //header row
            for(int i = 0; i < Record.size(); i++) {
                for(int j = 0; j < 16; j++) {
                    test.write(String.valueOf(Record.get(i).get(j)) + ",");
                }
                test.write("\n");
            }
            test.close();
        }
        catch (IOException e) {
            System.out.println("an error occurred.");
        }
        
        
        
        //calculating instability
        int num = 0;
        double momentum_sum = 0.0;
        double init_momentum = Record.get(0).get(15);
        double stop_time = Record.get(Record.size() - 1).get(0);
        while((stop_time - Record.get(Record.size() - 1 - num).get(0)) < period) {
            momentum_sum += Record.get(Record.size() - 1 - num).get(15);
            num++;
        }
        double fin_momentum = momentum_sum / num;
        double instability = Math.abs(fin_momentum - init_momentum)/(Planet.mass * stop_time);
        
        double insta_err = 2*Math.pow(10,-12)*Math.pow(a, -1.782);
        //We use this error in instability when the stars' positions are "forced."
        //A different method is needed, the simulation parameters are changed.
        
        
        
        //writing this run to the master csv file
        try {
            //C:\Users\grimm\OneDrive\Documents\COLLEGE CLASSES\Senior Research\senior_research_code\
            FileWriter master = new FileWriter("master.csv", true);
            master.write(String.valueOf(Planet.mass) + ",");
            master.write(String.valueOf(periapsis) + ",");
            master.write(String.valueOf(eccentricity) + ",");
            master.write(String.valueOf(a) + ",");
            master.write(String.valueOf(period) + ",");
            master.write(String.valueOf(instability) + ",");
            master.write(String.valueOf(insta_err) + ",");
            
            master.write("\n");
            master.close();
        }
        catch (IOException e) {
            System.out.println("an error occurred.");
        }
       
    } // end of for loop 1
    } // end of for loop 2
    } // end of main ----------------------------------------------------------------------------------------------
    
    
    
    
    
    
    
    // Functions
    
    //Dynamic RK algorithms-----------------------------------------------------------------------------------------
    public static void dynamic_RK4(body b_1, body b_2, body b_3, double h) {
        b_1.RK4_jk1(b_2, b_3, h);
        b_2.RK4_jk1(b_1, b_3, h);
        b_3.RK4_jk1(b_1, b_2, h);
        
        b_1.RK4_jk2(b_2, b_3, h);
        b_2.RK4_jk2(b_1, b_3, h);
        b_3.RK4_jk2(b_1, b_2, h);
        
        b_1.RK4_jk3(b_2, b_3, h);
        b_2.RK4_jk3(b_1, b_3, h);
        b_3.RK4_jk3(b_1, b_2, h);
        
        b_1.RK4_jk4(b_2, b_3, h);
        b_2.RK4_jk4(b_1, b_3, h);
        b_3.RK4_jk4(b_1, b_2, h);
        
        b_1.set_RK4(h);
        b_2.set_RK4(h);
        b_3.set_RK4(h);
    }
    
    public static void dynamic_RK45(body b_1, body b_2, body b_3, double h, double tolerance) {
        b_1.RK45_jk1(b_2, b_3, h);
        b_2.RK45_jk1(b_1, b_3, h);
        b_3.RK45_jk1(b_1, b_2, h);
        
        b_1.RK45_jk2(b_2, b_3, h);
        b_2.RK45_jk2(b_1, b_3, h);
        b_3.RK45_jk2(b_1, b_2, h);
        
        b_1.RK45_jk3(b_2, b_3, h);
        b_2.RK45_jk3(b_1, b_3, h);
        b_3.RK45_jk3(b_1, b_2, h);
        
        b_1.RK45_jk4(b_2, b_3, h);
        b_2.RK45_jk4(b_1, b_3, h);
        b_3.RK45_jk4(b_1, b_2, h);
        
        b_1.RK45_jk5(b_2, b_3, h);
        b_2.RK45_jk5(b_1, b_3, h);
        b_3.RK45_jk5(b_1, b_2, h);
        
        b_1.RK45_jk6(b_2, b_3, h);
        b_2.RK45_jk6(b_1, b_3, h);
        b_3.RK45_jk6(b_1, b_2, h);
        
        b_1.set_sol4(h); b_1.set_sol5(h);
        b_2.set_sol4(h); b_2.set_sol5(h);
        b_3.set_sol4(h); b_3.set_sol5(h);
    }
    
    
    
    
    // other functions --------------------------------------------------------------------------------------------
    //adapts the new time step given a tolerance.
    public static double adapt_Step(body b, double h, double tolerance) { 
        double diff = vec_Magnitude(vec_Diff(b.sol4_position, b.sol5_position));
        double h_prime = h;
        if(!(diff == 0)) {
            h_prime = h*(Math.pow((tolerance / diff), 0.2));
        }
        return h_prime;
    }
    
    public static double total_Energy(body b_1, body b_2, body b_3) {
        double energy = 0;
        energy += -vec_Magnitude(vec_Scalmult(b_1.Gravity(b_2, b_1.current_position), b_1.mass*vec_Magnitude(vec_Diff(b_1.current_position,b_2.current_position))));
        energy += -vec_Magnitude(vec_Scalmult(b_1.Gravity(b_3, b_1.current_position), b_1.mass*vec_Magnitude(vec_Diff(b_1.current_position,b_3.current_position))));
        energy += -vec_Magnitude(vec_Scalmult(b_2.Gravity(b_3, b_2.current_position), b_2.mass*vec_Magnitude(vec_Diff(b_2.current_position,b_3.current_position))));
        
        energy += 0.5*b_1.mass*Math.pow(vec_Magnitude(b_1.current_velocity), 2);
        energy += 0.5*b_2.mass*Math.pow(vec_Magnitude(b_2.current_velocity), 2);
        energy += 0.5*b_3.mass*Math.pow(vec_Magnitude(b_3.current_velocity), 2);
        return energy;
    }
    
    public static ArrayList<Double> center_Of_Mass(body b_1, body b_2, body b_3) {
        ArrayList<Double> x_1 = vec_Scalmult(b_1.current_position, b_1.mass);
        ArrayList<Double> x_2 = vec_Scalmult(b_2.current_position, b_2.mass);
        ArrayList<Double> x_3 = vec_Scalmult(b_3.current_position, b_3.mass);
        double tot_m = b_1.mass + b_2.mass + b_3.mass;
        return vec_Scalmult(vec_Add(x_1, vec_Add(x_2, x_3)), 1/tot_m);
    }
    
    public static ArrayList<Double> center_Of_Mass_Vel(body b_1, body b_2, body b_3) {
        ArrayList<Double> v_1 = vec_Scalmult(b_1.current_velocity, b_1.mass);
        ArrayList<Double> v_2 = vec_Scalmult(b_2.current_velocity, b_2.mass);
        ArrayList<Double> v_3 = vec_Scalmult(b_3.current_velocity, b_3.mass);
        double tot_m = b_1.mass + b_2.mass + b_3.mass;
        return vec_Scalmult(vec_Add(v_1, vec_Add(v_2, v_3)), 1/tot_m);
    }
    
    //sets the center of mass to 0,0
    public static void com_Correction(body b_1, body b_2, body b_3) {
        ArrayList<Double> com = center_Of_Mass(b_1, b_2, b_3);
        b_1.current_position = vec_Diff(b_1.current_position, com);
        b_2.current_position = vec_Diff(b_2.current_position, com);
        b_3.current_position = vec_Diff(b_3.current_position, com);
    }
    
    //sets the center of mass velocity to 0,0
    public static void com_Vel_Correction(body b_1, body b_2, body b_3) {
        ArrayList<Double> com_vel = center_Of_Mass_Vel(b_1, b_2, b_3);
        b_1.current_velocity = vec_Diff(b_1.current_velocity, com_vel);
        b_2.current_velocity = vec_Diff(b_2.current_velocity, com_vel);
        b_3.current_velocity = vec_Diff(b_3.current_velocity, com_vel);
    }
    
    //constrains energy to E0 by scaling all bodies by the same factor.
    //WARNING: if there is a sign change in total energy, this function will produce inaccurate results. The simulation will still run.
    public static void constrain_Energy(body b_1, body b_2, body b_3, double E0) {
        double A = Math.abs(E0/total_Energy(b_1, b_2, b_3));
        b_1.current_velocity = vec_Scalmult(b_1.current_velocity, Math.sqrt(A));
        b_2.current_velocity= vec_Scalmult(b_2.current_velocity, Math.sqrt(A));
        b_3.current_velocity = vec_Scalmult(b_3.current_velocity, Math.sqrt(A));
        b_1.current_position = vec_Scalmult(b_1.current_position, 1/(A));
        b_2.current_position = vec_Scalmult(b_2.current_position, 1/(A));
        b_3.current_position = vec_Scalmult(b_3.current_position, 1/(A));
    }
    
    
    
    
    
    // Vector Operations --------------------------------------------------------------------------------------------
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
