package LifeCycleProblem;

import Jama.Matrix;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class Testing {

    double mortgage_rate(double H_0, double r_0, double mu_H, double mu_r,
                         double sigma_H, double sigma_r, double rho,
                         double T, int N, double OLTV){
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

        // Create an empty schedule for balance and add the initial balance (principal) to the schedule
        ArrayList<Double> balance = new ArrayList<>();
        balance.add(0, H_0 * OLTV);

        /* The critical steps to calculate the the mortgage rate including the main lambda function that creates
         * a mapping from a given mortgage rate (a guess) -> value of the mortgage at the origin node (given by
         * V(N, N)) given all the other parameters including probabilities of default and prepayment as a function
         * of FICO and LTV (origin and instantaneous).
         * */
        UnivariateFunction function = annual_mortgage_rate -> {

            // first convert the annual mortgage rate to monthly
            double mortgage_rate = annual_mortgage_rate / 12;

            // now calculate the monthly mortgage payment, x.
            double x = (mortgage_rate * balance.get(0)) / (1 - pow((1 + mortgage_rate), -N));
            System.out.println(x);

            // now calculate the complete balance schedule
            for (int i = 1; i <= N; i++) {
                balance.add(i, (1 + mortgage_rate) * balance.get(i - 1) - x);
            }


            // Backward induct
            for (int i = N - 1; i >= 0; i--) {
                for (int j = -i + N; j <= i + N; j += 2) {
                    for (int k = -i + N; k <= i + N; k += 2) {
                        V.set(j, k, exp(-rt.get(k) * dt) * (x + (pdd * V.get(j - 1, k - 1)
                                + pud * V.get(j + 1, k - 1) + pdu * V.get(j - 1, k + 1)
                                + puu * V.get(j + 1, k + 1))));
                        // Use this for rational behavior
                        V.set(j, k, min(V.get(j, k), min(balance.get(i), Ht.get(j))));
                    }
                }
            }

            /* calculate and return the difference between the principal of the loan and it's value as
             * a function of the mortgage rate (V(m | FICO, OLTV))
             * */
            System.out.println(balance.get(0) + " , " +  V.get(N, N));
            return balance.get(0) - V.get(N, N);
        };

        double error = function.value(0.2);
        // Instantiate the root finding algorithm from Apache Commons
        // UnivariateSolver solver = new BisectionSolver(absolute_accuracy);    // Bisection algo
        //UnivariateSolver solver = new PegasusSolver(absolute_accuracy);         // Pegasus algo
        // find and return the root (i.e. the (annual) mortgage rate such that value of the loan equals the principal)
        //return max(solver.solve(number_of_iterations, function, min_m, max_m, start_m), r_0);
        return error;
        // end of mortgage rate finding process; i.e, a single point on the surface has been found at this stage.
    }


    public static void main(String[] args){
        long startTime = System.nanoTime();

        Grids grids = new Grids();      // call the Grids class; it contains the constructor for the matrices
        String path = "/home/ahyan/Dropbox/Housing Market and Leverage Cycle/Code/src/LifeCycleProblem/Data/";
        // Params for the entire credit surface
        double[] parameters = grids.Parameters(path + "ParametersForCS.csv");
        int H_0 = (int) parameters[0];   // house price index at origin
        double r_0 = parameters[1];      // interest rate (discount rate) at origin
        double mu_H = parameters[2];     // drift for house price index
        double mu_r = parameters[3];     // drift for the interest rate
        double sigma_H = parameters[4];  // vol for house price index
        double sigma_r = parameters[5];  // vol for interest rate
        double rho = parameters[6];      // corr. between house price index and interest rate processes
        double T = parameters[7];        // life of a mortgage (years)
        int N = (int) parameters[8];     // life of a mortgage (months)

        // Params for the root finding algorithm (bisections)
        double min_m = r_0 - 5e-3;          // minimum possible value for the mortgage rate, m
        double max_m = r_0 + 0.2;           // maximum possible value for the mortgage rate, m
        double start_m = r_0;               // starting point for the mortgage rate search (always start at the risk-free rate)
        double absolute_accuracy = parameters[9];            // degree of accuracy desired for the root finding process (smaller number -> more run time)
        int number_of_iterations = (int) parameters[10];     // max number of iterations allowed (it may not use all iterations, this is a ceiling)

        // Params for the probabilities of prepayment and default matrices
        int numberOfFicoBins = (int) parameters[11];       // number of FICO bins for credit surface
        int numberOfLTVBins = (int) parameters[12];        // number of ltv bins for credit surface

        // completing the path for probabilities
        String defaultRatesPath = path + "DefaultRates.csv";
        String prepayRatesPath  = path + "PrepayRates.csv";

        Matrix DefaultRates = grids.DefaultRates(numberOfFicoBins, numberOfLTVBins, defaultRatesPath);
        Matrix PrepayRates = grids.PrepayRates(numberOfFicoBins, numberOfLTVBins, prepayRatesPath);

        // FICOs and OLTVs for which to calculate the mortgage rate
        double[] FICO = {500, 550, 600, 650, 700, 750, 800};
        double[] OLTV = {0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0};

        // 2D array to capture the mortgage rate for aforementioned FICOs and OLTVs
        double[][] credit_surface = new double[FICO.length][OLTV.length];

        // create the Cartesian set of OLTV x FICO
        Integer[][] matrix1 = new Integer[2][];
        matrix1[0] = ArrayUtils.toObject(IntStream.range(0, OLTV.length).toArray());
        matrix1[1] = ArrayUtils.toObject(IntStream.range(0, FICO.length).toArray());

        CartesianSet<Integer> stateSpace1 = new CartesianSet<>(matrix1);
        int Z = (int) stateSpace1.getCount();

        Testing testing = new Testing();
        double error = testing.mortgage_rate(H_0, r_0, mu_H, mu_r, sigma_H, sigma_r, rho, T, N, 0.8);


    }
}
