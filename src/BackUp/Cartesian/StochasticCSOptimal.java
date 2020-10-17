package BackUp;/* This file develops a two factor credit surface with GBMs for house price and interest rates.
This structure may be considered a prototype:
in the main body, it has optimal behaviour (min function).
I have removed the importing of files and other related variables so that when the model is iterated
it doesn't have to bring in the files again and again.
*/

import Jama.Matrix;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.javatuples.Triplet;

import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

import static java.lang.Math.*;


public class StochasticCSOptimal {

    public Triplet<Double, Double, Integer> CS_point(double H_0, double r_0, double mu_H, double mu_r,
                                                     double sigma_H, double sigma_r, double rho,
                                                     double T, int N, double OLTV, int m_index_guess){



        // set up mortgage grids
        int numberOfMortgageGridPoints = 2000-1;
        double stepOfMortgageGrid = 0.0001;

        // Set the stochastic process terms
        double dt = T / N;
        double dx_H = sigma_H * sqrt(dt);
        double dx_r = sigma_r * sqrt(dt);

        //Compute the probabilities
        double puu = ((dx_H * dx_r + (dx_r * mu_H + dx_H * mu_r + rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));
        double pud = ((dx_H * dx_r + (dx_r * mu_H - dx_H * mu_r - rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));
        double pdu = ((dx_H * dx_r - (dx_r * mu_H - dx_H * mu_r + rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));
        double pdd = ((dx_H * dx_r - (dx_r * mu_H + dx_H * mu_r - rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));

        ArrayList<Double> Ht = new ArrayList<>(Collections.nCopies(2 * N + 1, 0d));
        ArrayList<Double> rt = new ArrayList<>(Collections.nCopies(2 * N + 1, 0d));
        Matrix V = new Matrix(2 * N + 1, 2 * N + 1);

        // Initialize asset prices at (worst) maturity
        Ht.set(0, H_0 * exp(-N * dx_H));
        rt.set(0, r_0 * exp(-N * dx_r));

        // Compute stock prices at each node
        for (int j = 1; j <= 2 * N; j++){
            Ht.set(j, Ht.get(j - 1) * exp(dx_H));
            rt.set(j, rt.get(j - 1) * exp(dx_r));
        }

        // Set mortgage rate grid
        ArrayList<Double> mortgageRateGrid = new ArrayList<>();
        mortgageRateGrid.add(0, r_0);
        for (int i = 1; i < numberOfMortgageGridPoints + 1; i++){mortgageRateGrid.add(i, mortgageRateGrid.get(i - 1) + stepOfMortgageGrid);}
        //System.out.println("initial grid point = " + interestRateGrid.get(0) + " terminal grid point = " + interestRateGrid.get(numberOfInterestGridPoints));



        double[][] V_print = V.getArray();
        //for (int i = 0; i < V.getRowDimension(); i++){System.out.println(Arrays.toString(V_print[i]));}

        // Set the mortgage terms
        ArrayList<Double> balance = new ArrayList<>();
        balance.add(0, H_0 * OLTV);

        double mortgage_rate = 0d;
        int m_index = 0;

        for(int m = m_index_guess; m < mortgageRateGrid.size(); m++) {
            mortgage_rate = mortgageRateGrid.get(m) / 12;
            double x = (mortgage_rate * balance.get(0)) / (1 - pow((1 + mortgage_rate), -N));

            for (int i = 1; i <= N; i++) {
                balance.add(i, (1 + mortgage_rate) * balance.get(i - 1) - x);
            }

            //System.out.println("x = " + x);
            //System.out.println("balance = " + balance + "\n");


            double ltv = 0d;
            int l = 0;

            // Backward induct
            for (int i = N - 1; i >= 0; i--) {
                for (int j = -i + N; j <= i + N; j += 2) {
                    for (int k = -i + N; k <= i + N; k += 2) {
                        V.set(j, k, exp(-rt.get(k) * dt) * (x + (pdd * V.get(j - 1, k - 1)
                                + pud * V.get(j + 1, k - 1) + pdu * V.get(j - 1, k + 1)
                                + puu * V.get(j + 1, k + 1))));
                        V.set(j, k, min(V.get(j, k), min(balance.get(i), Ht.get(j))));
                        //V.set(j, k, prob_prepay * balance.get(i) + prob_default * balance.get(i) * recoveryFraction + (1 - prob_prepay - prob_default) * V.get(j, k));
                    }
                }
            }

            //for (int i = 0; i < V.getRowDimension(); i++){System.out.println(Arrays.toString(V_print[i]));}

            DecimalFormat decimalFormat = new DecimalFormat("0.00");
            //System.out.println("Error = " + abs(V.get(N, N) - balance.get(0)));


            // Compute the error and if it is less than tolerance then end the run
            if (abs(balance.get(0) - V.get(N, N)) <= 0.025){
                //System.out.println("OLTV = " + OLTV + ", m = " + decimalFormat.format(mortgage_rate * 12 * 100));
                m_index = m;
                break;
            }

        }

        Triplet<Double, Double, Integer> CS_point = new Triplet<>(OLTV, mortgage_rate * 12 * 100, m_index);


        return CS_point;

    }



    public static void main(String args[]){

        long startTime = System.nanoTime();

        DecimalFormat decimalFormat = new DecimalFormat("0.00");

        // the OLTV points for which I want to calculate the mortgage rate
        double[] OLTV = {0.3, 0.6, 0.7, 0.75, 0.8, 0.825, 0.875, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0};

        // list for capturing the mortgage rates for the credit surface
        double[] M = new double[OLTV.length];

        // inital guess for mortgage rate index
        int m_index_guess = 0;

        // loop through CS_point for each OLTV to find the corresponding mortgage rate and it's index
        for (int i = 0; i < OLTV.length; i  ++){
            Triplet<Double, Double, Integer> CS_point = new StochasticCSOptimal().CS_point(125, 0.03, 0, 0, 0.0623, 0.0031,
                    0, 30, 360, OLTV[i], m_index_guess);

            System.out.println("OLTV = " + CS_point.getValue0() + " , mortgage_rate = " + CS_point.getValue1() + " , m_index = " + CS_point.getValue2());

            M[i] = CS_point.getValue1();

            m_index_guess = CS_point.getValue2();
        }


        double[] OLTV_interpolated = {0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
                0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1.0};
        double[] M_interpolated = new double[OLTV_interpolated.length];

        UnivariateInterpolator univariateInterpolator = new LinearInterpolator();

        UnivariateFunction univariateFunction = univariateInterpolator.interpolate(OLTV, M);

        for (int i = 0; i < OLTV_interpolated.length; i++){
            M_interpolated[i] = univariateFunction.value(OLTV_interpolated[i]);
        }

        // print the Credit Surface to .csv file
        try(FileWriter fileWriter = new FileWriter("C:\\Users\\ahyan\\Dropbox\\Housing Market and Leverage Cycle\\Code\\StochCSMin.csv")){
            // Write header first
            fileWriter.append("OLTV");
            fileWriter.append(",");
            fileWriter.append("m (%)");
            fileWriter.append("\n");
            // Loop through to write the rest of the CS
            for (int i = 0; i < OLTV_interpolated.length; i++){
                fileWriter.append(String.valueOf(OLTV_interpolated[i]));
                fileWriter.append(",");
                fileWriter.append(String.valueOf(M_interpolated[i]));
                fileWriter.append("\n");
            }
            fileWriter.close();
        }catch (IOException ioe){ioe.printStackTrace();}



        long endTime = System.nanoTime();

        long duration = (endTime - startTime);

        System.out.println(Math.round(duration * 1e-9) + " seconds");
    }
}
