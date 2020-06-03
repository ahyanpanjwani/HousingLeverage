/* This file develops a two factor credit surface with GBMs for house price and interest rates.
This structure may be considered a prototype:
in the main body, it has realistic behavior in that we use probs of
prepay and default or rational behavior (min fn) depending on what is commented in/out.
*/

import Jama.Matrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.Math.*;

public class StochasticCSRealistic {

    public double mortgage_rate(double H_0, double r_0, double mu_H, double mu_r, double sigma_H, double sigma_r, double rho,
                        double T, int N, double OLTV, int FICO_index){

        // Import probabilities

        int numberOfFicoBins = 7;
        int numberOfLTVBins = 12;
        int numberOfMortgageGridPoints = 2000-1;
        double stepOfMortgageGrid = 0.0001;

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
        double[][] defaultArrays = DefaultRates.getArray();
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
        double[][] prepayArrays = PrepayRates.getArray();
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
                        V.set(j, k, min(V.get(j, k), min(balance.get(i), Ht.get(j))));
                        //V.set(j, k, prob_prepay * balance.get(i) + prob_default * balance.get(i) * recoveryFraction + (1 - prob_prepay - prob_default) * V.get(j, k));
                    }
                }
            }

            //for (int i = 0; i < V.getRowDimension(); i++){System.out.println(Arrays.toString(V_print[i]));}

            DecimalFormat decimalFormat = new DecimalFormat("0.00");
            System.out.println("Error = " + abs(V.get(N, N) - balance.get(0)));


            if (abs(balance.get(0) - V.get(N, N)) <= 0.025){
                System.out.println("Principal = " + balance.get(0));
                System.out.println("Value = " + V.get(N, N));
                System.out.println("Error = " + abs(V.get(N, N) - balance.get(0)));
                System.out.println("mortgage rate = " + decimalFormat.format(mortgage_rate * 12 * 100) + "%");
                System.out.println("m_index = " + m);
                break;
            }

        }

        return mortgage_rate * 12;
    }



    public static void main(String args[]){
        double mortgage_rate = new StochasticCSRealistic().mortgage_rate(125, 0.03, 0.0823, -0.0001, 0.0623, 0.0031,
                0, 30, 360, 1.0, 0);
    }
}
