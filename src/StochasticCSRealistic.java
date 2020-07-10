/* This file develops a two factor credit surface with GBMs for house price and interest rates.
This structure may be considered a prototype:
in the main body, it has realistic behavior in that we use probs of
prepay and default or rational behavior (min fn) depending on what is commented in/out.
*/

import Jama.Matrix;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.javatuples.Triplet;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class StochasticCSRealistic {

    public Triplet<Double, Double, Integer> CS_point(double H_0, double r_0, double mu_H, double mu_r,
                                                     double sigma_H, double sigma_r, double rho,
                                                     double T, int N, double OLTV, int FICO_index,
                                                     int m_index_guess, Matrix DefaultRates, Matrix PrepayRates){

        //---------SET UP MORTGAGE GRIDS-----------//
        int numberOfMortgageGridPoints = 1000-1;
        double stepOfMortgageGrid = 0.0002;

        // Set the stochastic process terms
        double dt = T / N;
        double dx_H = sigma_H * sqrt(dt);
        double dx_r = sigma_r * sqrt(dt);

        //Compute the probabilities
        double puu = ((dx_H * dx_r + (dx_r * mu_H + dx_H * mu_r + rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));
        double pud = ((dx_H * dx_r + (dx_r * mu_H - dx_H * mu_r - rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));
        double pdu = ((dx_H * dx_r - (dx_r * mu_H - dx_H * mu_r + rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));
        double pdd = ((dx_H * dx_r - (dx_r * mu_H + dx_H * mu_r - rho * sigma_H * sigma_r) * dt) / (4 * dx_H * dx_r));

        // Create grids to store trees
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
        //System.out.println("initial grid point = " + mortgageRateGrid.get(0) + " terminal grid point = " + mortgageRateGrid.get(numberOfMortgageGridPoints) + " len:" + mortgageRateGrid.size());
        //System.out.println(mortgageRateGrid);



        double[][] V_print = V.getArray();
        //for (int i = 0; i < V.getRowDimension(); i++){System.out.println(Arrays.toString(V_print[i]));}

        // Set the mortgage terms
        ArrayList<Double> balance = new ArrayList<>();
        balance.add(0, H_0 * OLTV);

        double mortgage_rate = 0d;
        int m_index = 0;

        for (int m = m_index_guess; m < mortgageRateGrid.size(); m++){
            mortgage_rate = mortgageRateGrid.get(m) / 12;
            double x = (mortgage_rate * balance.get(0)) / (1 - pow((1 + mortgage_rate), -N));

            for (int i = 1; i <= N; i++) {
                balance.add(i, (1 + mortgage_rate) * balance.get(i - 1) - x);
            }

            double ltv;
            int f = FICO_index;
            int l = 0;
            double prob_prepay;
            double prob_default;

            // Backward induct
            for (int i = N - 1; i >= 0; i--) {
                for (int j = -i + N; j <= i + N; j += 2) {
                    for (int k = -i + N; k <= i + N; k += 2) {
                        ltv = balance.get(i) / Ht.get(j);
                        double recoveryFraction = Math.min(1, 0.5 / ltv);
                        if (ltv <= 0.3) {
                            l = 0;
                        } else if (ltv > 0.3 && ltv <= 0.4) {
                            l = 1;
                        } else if (ltv > 0.4 && ltv <= 0.5) {
                            l = 2;
                        } else if (ltv > 0.5 && ltv <= 0.6) {
                            l = 3;
                        } else if (ltv > 0.6 && ltv <= 0.7) {
                            l = 4;
                        } else if (ltv > 0.7 && ltv <= 0.8) {
                            l = 5;
                        } else if (ltv > 0.8 && ltv <= 0.9) {
                            l = 6;
                        } else if (ltv > 0.9 && ltv <= 1.0) {
                            l = 7;
                        } else if (ltv > 1.0 && ltv <= 1.1) {
                            l = 8;
                        } else if (ltv > 1.1 && ltv <= 1.2) {
                            l = 9;
                        } else if (ltv > 1.2 && ltv <= 1.3) {
                            l = 10;
                        } else if (ltv > 1.3) {
                            l = 11;
                        }
                        prob_prepay = PrepayRates.get(f, l);
                        prob_default = DefaultRates.get(f, l);
                        V.set(j, k, exp(-rt.get(k) * dt) * (x + (pdd * V.get(j - 1, k - 1)
                                + pud * V.get(j + 1, k - 1) + pdu * V.get(j - 1, k + 1)
                                + puu * V.get(j + 1, k + 1))));
                        //V.set(j, k, min(V.get(j, k), min(balance.get(i), Ht.get(j))));
                        V.set(j, k, prob_prepay * balance.get(i) + prob_default * balance.get(i) * recoveryFraction + (1 - prob_prepay - prob_default) * V.get(j, k));
                    }
                }
            }

            // Compute the error and if it is less than tolerance then end the run
            if (abs(balance.get(0) - V.get(N, N)) <= 0.5){
                m_index = m;
                break;
            }
        }

        Triplet<Double, Double, Integer> CS_point = new Triplet<>(OLTV, mortgage_rate * 12 * 100, m_index);
        return CS_point;
    }



    public static void main(String args[]){

        long startTime = System.nanoTime();

        int numberOfFicoBins = 7;
        int numberOfLTVBins = 12;
        String defaultRatesPath = "DataWork/LeverageCycle/DefaultRates.csv";
        String prepayRatesPath = "DataWork/LeverageCycle/PrepayRates.csv";

        Grids grids = new Grids();

        Matrix DefaultRates = grids.DefaultRates(numberOfFicoBins, numberOfLTVBins, defaultRatesPath);
        Matrix PrepayRates = grids.PrepayRates(numberOfFicoBins, numberOfLTVBins, prepayRatesPath);

        StochasticCSRealistic stochasticCSRealistic = new StochasticCSRealistic();

        double[] FICO = {500, 550, 600, 650, 700, 750, 800};
        double[] OLTV = {0.3, 0.6, 0.7, 0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0};

        double[][] credit_surface = new double[FICO.length][OLTV.length];

        AtomicInteger m_index_guess = new AtomicInteger();

        int H_0 = 60;

        IntStream.range(0, FICO.length).parallel().forEach(z ->{
            System.out.println("FICO: " + FICO[z]);
            for (int ioltv = 0; ioltv < OLTV.length; ioltv++){
                Triplet<Double, Double, Integer> CS_point = stochasticCSRealistic.CS_point(H_0,
                        0.03, 0, 0, 0.0623, 0.0031, 0, 30, 360,
                        OLTV[ioltv], z, m_index_guess.get(), DefaultRates, PrepayRates);

                credit_surface[z][ioltv] = CS_point.getValue1();

                m_index_guess.set(CS_point.getValue2());
            }

        });

        double[] OLTV_interpolated = {0.3, 0.6, 0.65, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.84, 0.86, 0.88,
                0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0};
        double[][] credit_surface_interpolated = new double[FICO.length][OLTV_interpolated.length];

        UnivariateInterpolator univariateInterpolator =  new LinearInterpolator();

        for (int f = 0; f < FICO.length; f++){

            UnivariateFunction univariateFunction = univariateInterpolator.interpolate(OLTV, credit_surface[f]);

            for (int l = 0; l < OLTV_interpolated.length; l++){
                credit_surface_interpolated[f][l] = univariateFunction.value(OLTV_interpolated[l]);
            }
        }

        try(FileWriter fileWriter = new FileWriter("DataWork/LeverageCycle/StochCSRealistic" + String.valueOf(H_0) + ".csv")){
            for (int f = 0; f < FICO.length; f++){
                for (int l = 0; l < OLTV_interpolated.length; l++){
                    fileWriter.append(String.valueOf(FICO[f]));
                    fileWriter.append(",");
                    fileWriter.append(String.valueOf(OLTV_interpolated[l]));
                    fileWriter.append(",");
                    fileWriter.append(String.valueOf(credit_surface_interpolated[f][l]));
                    fileWriter.append("\n");
                }
            }
            fileWriter.close();
        }catch (IOException ioe){ioe.printStackTrace();}


        long endTime = System.nanoTime();

        long duration = (endTime - startTime);

        System.out.println(round(duration * 1e-9) / 60 + " mins");
    }
}
