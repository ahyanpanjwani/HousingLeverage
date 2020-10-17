package Sparse;

import org.apache.commons.math3.distribution.NormalDistribution;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class Tauchen {

    // productivity grid, discretized a la Tauchen (1986)
    double[] productivityGrid(int ne, double sigma_eps, double lambda_eps, double m){
        double[] egrid = new double[ne];
        double sigma_y = sqrt(pow(sigma_eps, 2) / (1 - pow(lambda_eps, 2)));
        double estep = 2 * sigma_y * m / (ne - 1);
        egrid[0] = - m * sigma_y;
        for (int i = 1; i < egrid.length; i++){
            egrid[i] = egrid[i - 1] + estep;
        }

        return egrid;
    }

    NormalDistribution norm = new NormalDistribution();

    double[] transitionVector(double currentProductivity, double[] futureProductivity, int ne, double sigmaEps, double lambdaEps){

        double[] pi = new double[ne];
        double mm = futureProductivity[1] - futureProductivity[0];
        for (int k = 0; k < ne; k++){
            if (k == 0){
                pi[k] = norm.cumulativeProbability((futureProductivity[k] - lambdaEps * currentProductivity + (mm / 2)) / sigmaEps);
            }else if (k == ne - 1){
                pi[k] = 1 - norm.cumulativeProbability((futureProductivity[k] - lambdaEps * currentProductivity - (mm / 2)) / sigmaEps);
            }else {
                pi[k] = norm.cumulativeProbability((futureProductivity[k] - lambdaEps * currentProductivity + (mm/2))/sigmaEps)
                        - norm.cumulativeProbability((futureProductivity[k] - lambdaEps * currentProductivity - (mm/2))/sigmaEps);
            }
        }

        return pi;
    }

}
