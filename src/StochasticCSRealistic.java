/* This file develops a two factor credit surface with GBMs for house price and interest rates.
This structure may be considered a prototype:
in the main body, it has realistic behavior in that we use probs of
prepay and default or rational behavior (min fn) depending on what is commented in/out.
*/

import Jama.Matrix;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class StochasticCSRealistic {

    public double[] mortgage_rates(double H_0, double r_0, double mu_H, double mu_r, double sigma_H, double sigma_r, double rho,
                                   double T, int N, int FICO_index){

        // Import probabilities

        int numberOfFicoBins = 7;
        int numberOfLTVBins = 12;
        int numberOfMortgageGridPoints = 1000-1;
        double stepOfMortgageGrid = 0.0002;

        String defaultRatesPath = "C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Data\\MonthlyDefaultData\\" +
                "DefaultRates" + "Jan-02" + ".csv";
        String prepayRatesPath = "C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Data\\MonthlyPrepayData\\" +
                "PrepayRates" + "Jan-02" + ".csv";
        /////////DEFAULT MATRIX////////////////
        List<String[]> rowList = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(defaultRatesPath))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList.add(lineItems);
            }
            br.close();
        }catch (Exception e){e.printStackTrace();}


        String[][] Default = new String[rowList.size()][];
        for (int i = 0; i < rowList.size(); i++){
            String[] row = rowList.get(i);
            Default[i] = row;
        }


        Matrix DefaultRates = new Matrix(numberOfFicoBins, numberOfLTVBins);
        for (int j = 0; j < DefaultRates.getRowDimension(); j++) {
            for (int i = 0; i < DefaultRates.getColumnDimension(); i++) {
                DefaultRates.set(j, i, Double.parseDouble(Default[j+1][i + 1]));
            }
        }

        //Loop through default matrix to print:
        //double[][] defaultArrays = DefaultRates.getArray();
        //for (int i = 0; i < DefaultRates.getRowDimension(); i++){ System.out.println(Arrays.toString(defaultArrays[i]));}
        //System.out.println("CUT");

        /////////PREPAY MATRIX////////////////
        List<String[]> rowList2 = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(prepayRatesPath))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList2.add(lineItems);
            }
            br.close();
        }catch (Exception e){}


        String[][] Prepay = new String[rowList2.size()][];
        for (int i = 0; i < rowList2.size(); i++){
            String[] row = rowList2.get(i);
            Prepay[i] = row;
        }



        Matrix PrepayRates = new Matrix(numberOfFicoBins, numberOfLTVBins);
        for (int j = 0; j < PrepayRates.getRowDimension(); j++) {
            for (int i = 0; i < PrepayRates.getColumnDimension(); i++) {
                PrepayRates.set(j, i, Double.parseDouble(Prepay[j+1][i + 1]));
            }
        }

        //Loop through prepay matrix to print:
        //double[][] prepayArrays = PrepayRates.getArray();
        //for (int i = 0; i < PrepayRates.getRowDimension(); i++){ System.out.println(Arrays.toString(prepayArrays[i]));}
        //System.out.println("CUT");


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
        //System.out.println("initial grid point = " + mortgageRateGrid.get(0) + " terminal grid point = " + mortgageRateGrid.get(numberOfMortgageGridPoints));
        //System.out.println(mortgageRateGrid);



        double[][] V_print = V.getArray();
        //for (int i = 0; i < V.getRowDimension(); i++){System.out.println(Arrays.toString(V_print[i]));}

        double mortgage_rate = 0d;
        int m_index_guess = 0;
        int m_index = 0;

        double[] OLTV = {0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.88, 0.9, 0.92, 0.95, 0.98, 1.0}; // placeholder
        double[] M = new double[OLTV.length];

        for (int o = 0; o < OLTV.length; o++) {

            // Set the mortgage terms
            double oltv = OLTV[o];
            ArrayList<Double> balance = new ArrayList<>();
            balance.add(0, H_0 * oltv);

            for (int m = m_index_guess; m < mortgageRateGrid.size(); m++) {
                mortgage_rate = mortgageRateGrid.get(m) / 12;
                double x = (mortgage_rate * balance.get(0)) / (1 - pow((1 + mortgage_rate), -N));

                for (int i = 1; i <= N; i++) {
                    balance.add(i, (1 + mortgage_rate) * balance.get(i - 1) - x);
                }
                //System.out.println("x = " + x);
                //System.out.println("balance = " + balance + "\n");


                double ltv = 0d;
                int f = FICO_index;
                int l = 0;
                double prob_prepay = 0;
                double prob_default = 0;

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

                //for (int i = 0; i < V.getRowDimension(); i++){System.out.println(Arrays.toString(V_print[i]));}

                DecimalFormat decimalFormat = new DecimalFormat("0.00");
                //System.out.println("Error = " + abs(V.get(N, N) - balance.get(0)));


                if (abs(balance.get(0) - V.get(N, N)) <= 0.5) {
                    //System.out.println("Principal = " + balance.get(0));
                    //System.out.println("Value = " + V.get(N, N));
                    //System.out.println("Error = " + abs(V.get(N, N) - balance.get(0)));
                    //System.out.println("mortgage rate = " + decimalFormat.format(mortgage_rate * 12 * 100) + "%");
                    //System.out.println("m_index = " + m);
                    m_index = m;
                    break;
                }

            }

            M[o] = mortgage_rate * 12 * 100;
            m_index_guess = m_index;
            //System.out.println("OLTV = " + oltv + ", mortgage rate = " + mortgage_rate * 12 * 100 + ", m_index_guess = " + m_index);

        }

        UnivariateInterpolator univariateInterpolator = new LinearInterpolator();

        UnivariateFunction univariateFunction = univariateInterpolator.interpolate(OLTV, M);

        double[] OLTV_interpolated = {0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.84, 0.86, 0.88,
                0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0};
        double[] M_interpolated = new double[OLTV_interpolated.length];

        for (int i = 0; i < OLTV_interpolated.length; i++){
            M_interpolated[i] = univariateFunction.value(OLTV_interpolated[i]);
        }

        return M_interpolated;
    }



    public static void main(String args[]){

        long startTime = System.nanoTime();

        double[] FICO = {500, 550, 600, 650, 700, 750, 800};
        double[] OLTV_interpolated = {0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.84, 0.86, 0.88,
                0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0};

        double[][] credit_surface = new double[FICO.length][OLTV_interpolated.length];


        StochasticCSRealistic stochasticCSRealistic = new StochasticCSRealistic();

        IntStream.range(0, FICO.length).parallel().forEach(z -> {
            System.out.println("FICO = " + FICO[z]);
            credit_surface[z] = stochasticCSRealistic.mortgage_rates(125, 0.03, 0.00, 0, 0.0623, 0.0031, 0,
                    30, 360, z);
        });



        try(FileWriter fileWriter = new FileWriter("C:\\Users\\ahyan\\Dropbox\\Housing Market and Leverage Cycle\\Code\\StochCSRealistic.csv")){
            for (int i = 0; i < FICO.length; i++){
                for (int j = 0; j < OLTV_interpolated.length; j++){
                    fileWriter.append(String.valueOf(FICO[i]));
                    fileWriter.append(",");
                    fileWriter.append(String.valueOf(OLTV_interpolated[j]));
                    fileWriter.append(",");
                    fileWriter.append(String.valueOf(credit_surface[i][j]));
                    fileWriter.append("\n");
                }
            }
            fileWriter.close();
        }catch (IOException ioe){ioe.printStackTrace();}




        long endTime = System.nanoTime();

        long duration = (endTime - startTime);

        System.out.println(Math.round(duration * 1e-9) + " seconds");
    }
}
