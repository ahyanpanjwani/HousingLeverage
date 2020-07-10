import Jama.Matrix;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

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
    double[][][][][][] V(int T, int na, int ne, int nq, int nCS, int nn){
        double[][][][][][] V = new double[T][na][ne][nq][nCS][nn];
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

}
