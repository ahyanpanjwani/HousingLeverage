package BackUp;

import Jama.Matrix;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static java.lang.Math.*;

/**
 * This file creates the necessary grids (for a, e, q, CS, n) so that I don't have to create them for each model iteration.
 */

public class Grids {

    // asset grid
    double[] agrid(int na, double amin, double amax){

        double[] agrid = new double[na];
        agrid[0] = amin;
        double astep = (amax - amin) / (na - 1);

        for (int i = 1; i < agrid.length; i++){
            agrid[i] = agrid[i - 1] + astep;
        }

        return agrid;
    }

    // productivity grid, discretized a la Tauchen (1986)
    double[] egrid(int ne, double sigma_eps, double lambda_eps, double m){
        double[] egrid = new double[ne];
        double sigma_y = sqrt(pow(sigma_eps, 2) / (1 - pow(lambda_eps, 2)));
        double estep = 2 * sigma_y * m / (ne - 1);
        egrid[0] = - m * sigma_y;
        for (int i = 1; i < egrid.length; i++){
            egrid[i] = egrid[i - 1] + estep;
        }

        return egrid;
    }

    // productivity transition matrix (for the Markov transition)
    NormalDistribution norm = new NormalDistribution();
    Matrix P(int ne, double sigma_eps, double lambda_eps, double m, double[] egrid){
        Matrix P = new Matrix(ne, ne);
        double mm = egrid[1] - egrid[0];
        for (int j = 0; j < ne; j++){
            for (int k = 0; k < ne; k++){
                if (k == 0){
                    P.set(j, k, norm.cumulativeProbability((egrid[k] - lambda_eps*egrid[j] + (mm/2))/sigma_eps));
                }else if (k == ne - 1){
                    P.set(j, k, 1 - norm.cumulativeProbability((egrid[k] - lambda_eps*egrid[j] - (mm/2))/sigma_eps));
                }else {
                    P.set(j, k, norm.cumulativeProbability((egrid[k] - lambda_eps*egrid[j] + (mm/2))/sigma_eps)
                            - norm.cumulativeProbability((egrid[k] - lambda_eps*egrid[j] - (mm/2))/sigma_eps));
                }
            }
        }
        return P;
    }

    // grid for quantity of housing, q
    double[] qgrid(int nq, double qmin, double qmax){

        double[] qgrid = new double[nq];
        qgrid[0] = qmin;
        double qstep = (qmax - qmin) / (nq - 1);

        for (int i = 1; i < qgrid.length; i++){
            qgrid[i] = qgrid[i - 1] + qstep;
        }

        return qgrid;
    }

    // grid for oltv (it calls the credit surface csv and makes a grid for the oltvs for fico=500)
    double[] lgrid(String filePathForCS, int nl){
        double[] lgrid = new double[nl];
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(filePathForCS));
            String row;
            int it = 0;
            while ((row = bufferedReader.readLine()) != null & it < nl) {
                String[] data = row.split(",");
                lgrid[it] = Double.parseDouble(data[1]);
                it++;
            }
        }catch (IOException e){e.printStackTrace();}
        return lgrid;
    }

    // grid for mortgage rates for corresponding oltvs from the credit surface
    double[] mgrid(String filePathForCS, int nm){
        double[] mgrid = new double[nm];
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(filePathForCS));
            String row;
            int it = 0;
            while ((row = bufferedReader.readLine()) != null & it < nm) {
                String[] data = row.split(",");
                mgrid[it] = Double.parseDouble(data[2]) / 100;
                it++;
            }
        }catch (IOException e){e.printStackTrace();}
        return mgrid;
    }

    // grid for age of the mortgage used to calculate a households equity for when they want to sell/refinance/bequest
    double[] ngrid(int nn, double nmin, double nmax){
        double[] ngrid = new double[nn];
        ngrid[0] = nmin;
        double nstep = (nmax - nmin) / (nn - 1);

        for (int i = 1; i < ngrid.length; i++){
            ngrid[i] = ngrid[i - 1] + nstep;
        }

        return ngrid;
    }

    // the value function hypercube
    double[][][][][][] V(int T, int na, int nq, int nCS, int nl, int ne){
        //double[][][][][][] V = new double[T][na][ne][nq][nCS][nl];        // for use with model v.G and earlier
        double[][][][][][] V = new double[T][na][nq][nCS][nl][ne];          // for use with model v.H and later
        return V;
    }

    double[][] V_0(int na, int ne){
        double[][] V_0 = new double[na][ne];
        return V_0;
    }

    double[][][][][][] Expected(int T, int na, int nq, int nCS, int N, int ne, double[][][][][][] V, int age, Matrix P){
        double[][][][][][] Expected = new double[T][na][nq][nCS][N][ne];
        double expected = 0;
        for (int ia = 0; ia < na; ia++){
            for (int iq = 0; iq < nq; iq++){
                for (int iCS = 0; iCS < nCS; iCS++){
                    for (int in = 0; in < N; in++){
                        for (int ie = 0; ie < ne; ie++){
                            for (int iep = 0; iep < ne; iep++){
                                expected = expected + P.get(ie, iep) * V[age + 1][ia][iep][iq][iCS][in];
                            }
                            Expected[age + 1][ia][iq][iCS][in][ie] = expected;
                            expected = 0;
                        }
                    }
                }
            }
        }
        return Expected;
    }

    double[][][][][][] ExpectedB(int T, int na, int nq, int nCS, int nl, int N, int ne, double[][][][][][] V, int age,
                                 Matrix P, double p_h, double mu_h, double sigma_h,
                                 double[] mgrid, double[] lgrid, double[][] ltv_schedule){

        double[][][][][][] Expected = new double[T][na][nq][nCS][nl][ne];

        LinearInterpolator linearInterpolator = new LinearInterpolator();
        PolynomialSplineFunction polynomialSplineFunction;

        // price expectations terms (discretizing a GBM as a two-point process)
        double delta_x = sqrt(pow(sigma_h, 2) + pow(mu_h, 2));      // step for house price
        double pi_u = 0.5 + (0.5 * mu_h / delta_x);               // prob. of house price going up
        double pi_d = 0.5 - (0.5 * mu_h / delta_x);               // prob. of house price going down
        double p_h_u = p_h * exp(delta_x);                      // house price at the up node
        double p_h_d = p_h * exp(-delta_x);                     // house price at the down node

        for (int ia = 0; ia < na; ia++){
            for (int iq = 0; iq < nq; iq++){
                for (int iCS = 0; iCS < nCS; iCS++){
                    double mortgage_rate = mgrid[iCS];
                    double origination_ltv = lgrid[iCS];
                    double mortgage_payment = mortgage_rate * origination_ltv * p_h / (1 - pow(1 + mortgage_rate, - N));
                    for (int il = 0; il < nl; il++){
                        double current_ltv = ltv_schedule[iCS][il];
                        double balance_p = (1 + mortgage_rate) * current_ltv * p_h - mortgage_payment;
                        double ltv_u = min(ltv_schedule[iCS][0] - 1e-6, max(balance_p / p_h_u, 0 + 1e-6));
                        double ltv_d = min(ltv_schedule[iCS][0] - 1e-6, max(balance_p / p_h_d, 0 + 1e-6));
                        List<Double> ltv_list = Arrays.stream(ltv_schedule[iCS]).boxed().collect(Collectors.toList());
                        Collections.sort(ltv_list);
                        int idx_u = Collections.binarySearch(ltv_list, ltv_u);
                        int idx_d = Collections.binarySearch(ltv_list, ltv_d);
                        double[] ltv_bounds_u = new double[] {ltv_schedule[iCS][ltv_list.size() - abs(idx_u) + 1], ltv_schedule[iCS][ltv_list.size() - abs(idx_u)]};
                        double[] ltv_bounds_d = new double[] {ltv_schedule[iCS][ltv_list.size() - abs(idx_d) + 1], ltv_schedule[iCS][ltv_list.size() - abs(idx_d)]};
                        double[] expectation_bounds_u = new double[2];
                        double[] expectation_bounds_d = new double[2];
                        double expected_u, expected_d;
                        for (int ie = 0; ie < ne; ie++){
                            for (int iep = 0; iep < ne; iep++){
                                expectation_bounds_u[0] = expectation_bounds_u[0] + P.get(ie, iep) * V[age + 1][ia][iq][iCS][ltv_list.size() - abs(idx_u) + 1][iep];
                                expectation_bounds_u[1] = expectation_bounds_u[1] + P.get(ie, iep) * V[age + 1][ia][iq][iCS][ltv_list.size() - abs(idx_u)][iep];
                                expectation_bounds_d[0] = expectation_bounds_d[0] + P.get(ie, iep) * V[age + 1][ia][iq][iCS][ltv_list.size() - abs(idx_d) + 1][iep];
                                expectation_bounds_d[1] = expectation_bounds_d[1] + P.get(ie, iep) * V[age + 1][ia][iq][iCS][ltv_list.size() - abs(idx_d)][iep];
                            }
                            polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds_u, expectation_bounds_u);
                            expected_u = polynomialSplineFunction.value(ltv_u);
                            polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds_d, expectation_bounds_d);
                            expected_d = polynomialSplineFunction.value(ltv_d);
                            Expected[age + 1][ia][iq][iCS][il][ie] = pi_u * expected_u + pi_d * expected_d;
                        }
                    }
                }
            }
        }
        return Expected;
    }

    Matrix DefaultRates(int numberOfFicoBins, int numberOfLTVBins, String defaultRatesPath){
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

        return DefaultRates;
    }

    Matrix PrepayRates(int numberOfFicoBins, int numberOfLTVBins, String prepayRatesPath){
        List<String[]> rowList2 = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(prepayRatesPath))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList2.add(lineItems);
            }
            br.close();
        }catch (Exception e){e.printStackTrace();}


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
        return PrepayRates;
    }

    double[][] LTV_Schedule(double[] mgrid, double[] lgrid, int N, int nCS){
        double[][] ltv_schedule = new double[nCS][4 + N + 1];
        for (int iCS = 0; iCS < nCS; iCS++){
            for (int in = 0; in < N; in++){
                double origination_ltv = lgrid[iCS];
                double mortgage_rate = mgrid[iCS];
                ltv_schedule[iCS][in + 4] = origination_ltv * (pow(1 + mortgage_rate, N) - pow(1 + mortgage_rate, in)) / (pow(1 + mortgage_rate, N) - 1);
            }
            for (int in = 3; in >= 0; in--){
                double step = 0.05;
                ltv_schedule[iCS][in] = ltv_schedule[iCS][in + 1] + step;
            }
        }
        return ltv_schedule;
    }

}
