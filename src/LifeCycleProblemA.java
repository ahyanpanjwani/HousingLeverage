/**
 * This file is deprecated. It is an old model with buy only but the technique is from before I came across
 * Jesus Villaverde's code. LifeCycleProblemB has the new methods to solve lifecycle problems.
 **/

import com.google.common.primitives.Doubles;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Math.*;

public class LifeCycleProblemA {


    public LifeCycleProblemA(){
        int number_LTV_bins = 26;
        int number_FICO_bins = 7;

        // First read in the credit surface produced by StochasticCSRealistic

        String csvFile = "Dropbox\\Housing Market and Leverage Cycle\\Code\\StochCSRealistic.csv";
        List<String[]> rowList = new ArrayList<String[]>();
        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] lineItems = line.split(",");
                rowList.add(lineItems);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }


        String[][] cs = new String[rowList.size()][];
        for (int i = 0; i < rowList.size(); i++) {
            String[] row = rowList.get(i);
            cs[i] = row;
        }


        double[][] credit_surface = new double[number_FICO_bins][number_LTV_bins];
        for (int i = 0; i < number_FICO_bins; i++){
            for (int j = 0; j < number_LTV_bins; j++){
                credit_surface[i][j] = Double.parseDouble(cs[j + i * number_LTV_bins][2]) / 100;
            }
        }

        // SETTING UP THE REST OF THE STRUCTURE/GRIDS

        // Parameters
        double alpha = 0.33;
        double beta = 0.99;

        // Time Horizon
        int N = 30;

        // Prices
        double y = 0.01; // fixed annual wage
        double r_f = 3; // risk free rate
        double mu_H = 0; // drift for the housing process
        double p_0 = 125; // price guess
        double p_prime = p_0 * exp(mu_H);
        double q_supply = 1; // fixed supply of housing, normalized to a unit

        // Grids
        // assets
        double[] A = {1e-6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; //11 points, will probably have to increase

        //credit surface
        double[] OLTV = {0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.84, 0.86, 0.88,
                0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0}; //26 points (of interpolation)

        double[] M = credit_surface[0]; //mortgage rates for FICO = 500

        // housing quantity
        double[] Q = {0.0001, 0.00015, 0.0002, 0.00025, 0.0003, 0.00035, 0.0004, 0.00045, 0.0005, 0.00055, 0.0006,
                0.00065, 0.0007, 0.00075, 0.0008, 0.00085, 0.0009, 0.00095, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035}; //24 points

        // rental rate (per unit of quantity of housing rented). This number needs to be multiplied with q to get the (annual)
        // rent amount a renter shall pay. It assumes that landlord has an LTV of 0.8
        double rental_rate = M[10] * OLTV[10] * p_0 / (1 - pow(1 + M[10], -N)) + (p_0 - p_prime);

        // Last period renters problem
        double[][][] U_inter = new double[A.length][Q.length][A.length];
        double U_min = -1e6;

        double a_current = A[0];
        double a_future = A[0];

        double[] V_T_N_inter = new double[Q.length];

        for (int i = 0; i < Q.length; i++){
            double q = Q[i];
            double c = (1 + r_f) * a_current + y - a_future - rental_rate * q;
            double bequest = a_future;
            if (c > 0){
                V_T_N_inter[i] = (1 - alpha) * log(c) + alpha * log(q) + beta * log(bequest);
            }else {
                V_T_N_inter[i] = -1e6;
            }
        }

        System.out.println((1 + r_f) * a_current + y - a_future - (rental_rate * 0.003));

        System.out.println(Arrays.toString(V_T_N_inter));
        System.out.println(Doubles.max(V_T_N_inter));
        System.out.println(Q[Doubles.indexOf(V_T_N_inter, Doubles.max(V_T_N_inter))]);
        System.out.println(rental_rate);









    }


    public static void main(String args[]){
        LifeCycleProblemA lifeCycleProblemA = new LifeCycleProblemA();
    }



}
