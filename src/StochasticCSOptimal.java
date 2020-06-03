/* This file develops a two factor credit surface with GBMs for house price and interest rates.
This structure may be considered a prototype:
in the main body, it has optimal behaviour (min function).
I have removed the importing of files and other related variables so that when the model is iterated
it doesn't have to bring in the files again and again.
*/

import Jama.Matrix;
import javafx.util.Pair;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import static java.lang.Math.*;
import java.util.stream.IntStream;


public class StochasticCSOptimal {

    public double mortgage_rate(double H_0, double r_0, double mu_H, double mu_r, double sigma_H, double sigma_r, double rho,
                                double T, int N, double OLTV){

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

        for(int m = 0; m < mortgageRateGrid.size(); m++) {
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
                System.out.println("OLTV = " + OLTV + ", m = " + decimalFormat.format(mortgage_rate * 12 * 100));
                break;
            }

        }

        return mortgage_rate * 12 * 100;
    }



    public static void main(String args[]){

        long startTime = System.nanoTime();

        DecimalFormat decimalFormat = new DecimalFormat("0.00");

        // the OLTV points for which I want to calculate the mortgage rate
        double[] OLTV = {0.3, 0.6, 0.7, 0.75, 0.8, 0.825, 0.875, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0};

        ArrayList<String> M =new ArrayList<>(Collections.nCopies(OLTV.length, " "));
        ArrayList<Double> OLTV_processed = new ArrayList<>(Collections.nCopies(OLTV.length, 0d));

        // run the code in parallel given certain inputs
        IntStream.range(0, OLTV.length).parallel().forEach(z -> {
            double mortgage_rate = new StochasticCSOptimal().mortgage_rate(125, 0.03, 0, 0, 0.0623, 0.0031,
                    0, 30, 360, OLTV[z]);

            OLTV_processed.set(z, OLTV[z]);
            M.set(z, decimalFormat.format(mortgage_rate));
        });

        //
        System.out.println(OLTV_processed);
        System.out.println(M);

        // print the Credit Surface to .csv file
        try(FileWriter fileWriter = new FileWriter("C:\\Users\\ahyan\\Dropbox\\Housing Market and Leverage Cycle\\Code\\StochasticCreditSurface.csv")){
            // Write header first
            fileWriter.append("OLTV");
            fileWriter.append(",");
            fileWriter.append("m (%)");
            fileWriter.append("\n");
            // Loop through to write the rest of the CS
            for (int i = 0; i < M.size(); i++){
                fileWriter.append(String.valueOf(OLTV_processed.get(i)));
                fileWriter.append(",");
                fileWriter.append(M.get(i));
                fileWriter.append("\n");
            }
            fileWriter.close();
        }catch (IOException ioe){ioe.printStackTrace();}

        long endTime = System.nanoTime();

        long duration = (endTime - startTime);

        System.out.println(Math.round(duration * 1e-9) + " seconds");
    }
}
