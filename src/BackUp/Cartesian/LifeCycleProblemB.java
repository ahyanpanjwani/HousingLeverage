package BackUp; /**
 * This file replicates the Jesus Villaverde model for a simple life cycle model with assets and productivity
 * (the latter is stochastic). It also uses the Cartesian set technique to build/index the state space
 * for easier computation in parallel. This model does not contain any elements of the housing
 * model. Those will be in LifeCycleProblemC.java.
 */

import Jama.Matrix;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;
import java.util.stream.IntStream;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class LifeCycleProblemB {

    double[] xgrid(int nx, double xmin, double xmax){

        double[] xgrid = new double[nx];
        xgrid[0] = xmin;
        double xstep = (xmax - xmin) / (nx - 1);

        for (int i = 1; i < xgrid.length; i++){
            xgrid[i] = xgrid[i - 1] + xstep;
        }

        return xgrid;
    }

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

    double[][][] V(int T, int nx, int ne){
        double[][][] V = new double[T][nx][ne];
        return V; 
    }


    double[] value(int ie, int ix, int ne, int nx, int T, int age, double[] egrid, double[] xgrid, Matrix P,
               double[][][] V, double sigma, double beta, double w, double r){


        double expected;
        double utility;
        double consumption;
        double optConsumption = 0.0;
        double VV = pow(-10.0, 5.0);

        for (int ixp = 0; ixp < nx; ixp++){
            expected = 0;
            if (age < T - 1){
                for (int iep = 0; iep < ne; iep++){
                    expected = expected + P.get(ie, iep) * V[age + 1][ixp][iep];
                }
            }

            consumption = (1 + r) * xgrid[ix] + w * egrid[ie] - xgrid[ixp];

            utility = pow(consumption, 1-sigma) / (1-sigma) + beta*expected;

            if (consumption <= 0){utility = pow(-10, 5);}
            if (utility >= VV){
                VV = utility;
                optConsumption = consumption;
            }
        }

        double[] pair = new double[2];
        pair[0] = VV;
        pair[1] = optConsumption;
        return pair;
    }




    public static void main(String args[]){

        long startTime = System.nanoTime();

        int ne = 15;
        int nx = 1500;
        double xmin = 0.1;
        double xmax = 4.0;
        double sigma_eps = 0.02058;
        double lambda_eps = 0.99;
        double m = 1.5;
        int T = 10;
        double sigma = 2;
        double beta = 0.97;
        double w = 5;
        double r = 0.07;

        LifeCycleProblemB lifeCycleProblemB = new LifeCycleProblemB();
        double[] xgrid = lifeCycleProblemB.xgrid(nx, xmin, xmax);
        double[] egrid = lifeCycleProblemB.egrid(ne, sigma_eps, lambda_eps, m);
        Matrix P = lifeCycleProblemB.P(ne, sigma_eps, lambda_eps, m, egrid);
        double[][][] V = lifeCycleProblemB.V(T, nx, ne);
        double[][][] C = lifeCycleProblemB.V(T, nx, ne);

        Integer[][] matrix = new Integer[2][];
        matrix[0] = ArrayUtils.toObject(IntStream.range(0, nx).toArray());
        matrix[1] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());

        CartesianSet<Integer> stateSpace = new CartesianSet<>(matrix);

        for (int age = T - 1; age >= 0; age--){
            int finalAge = age;
            IntStream.range(0, nx * ne).parallel().forEach(z -> {
                V[finalAge][stateSpace.get(z).get(0)][stateSpace.get(z).get(1)] = lifeCycleProblemB.value(stateSpace.get(z).get(1),
                        stateSpace.get(z).get(0), ne, nx, T, finalAge, egrid, xgrid, P, V, sigma, beta, w, r)[0];
                C[finalAge][stateSpace.get(z).get(0)][stateSpace.get(z).get(1)] = lifeCycleProblemB.value(stateSpace.get(z).get(1),
                        stateSpace.get(z).get(0), ne, nx, T, finalAge, egrid, xgrid, P, V, sigma, beta, w, r)[1];
            });


            System.out.println("age = " + age);
        }
        System.out.println(Arrays.toString(C[0][800]));


        long endTime = System.nanoTime();
        System.out.println("run time = " + (endTime - startTime) * 1e-9);

    }

}
