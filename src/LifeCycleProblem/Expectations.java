package LifeCycleProblem;

import Jama.Matrix;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static java.lang.Math.*;
import static java.lang.Math.exp;

public class Expectations {

    public double[][][][][][][] ExpectationsWithFico(int T, int na, int nq, int nOLTV, int nLTV, int nFICO, int ne, int N,
                                                     double[][][][][][][] V, int age, Matrix P_E, Matrix P_F,
                                                     double p_h, double mu_h, double sigma_h, double[] lgrid,
                                                     double[][] creditSurface, double[][] ltvSchedule){
        double[][][][][][][] Expectation = new double[T][na][nq][nOLTV][nLTV][nFICO][ne];

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
                for (int iF = 0; iF < nFICO; iF++){
                    for (int iOlTV = 0; iOlTV < nOLTV; iOlTV++){
                        double mortgageRate = creditSurface[iF][iOlTV];
                        double originationLTV = lgrid[iOlTV];
                        double mortgagePayment = mortgageRate * originationLTV * p_h / (1 - pow(1 + mortgageRate, - N));
                        for (int iLTV = 0; iLTV < nLTV; iLTV++){
                            double currentLTV = ltvSchedule[iOlTV][iLTV];
                            double balancePrime = (1 + mortgageRate) * currentLTV * p_h - mortgagePayment;
                            double ltv_u = min(ltvSchedule[iOlTV][0] - 1e-6, max(balancePrime / p_h_u, 0 + 1e-6));
                            double ltv_d = min(ltvSchedule[iOlTV][0] - 1e-6, max(balancePrime / p_h_d, 0 + 1e-6));
                            List<Double> ltv_list = Arrays.stream(ltvSchedule[iOlTV]).boxed().collect(Collectors.toList());
                            Collections.sort(ltv_list);
                            int idx_u = Collections.binarySearch(ltv_list, ltv_u);
                            int idx_d = Collections.binarySearch(ltv_list, ltv_d);
                            double[] ltv_bounds_u = new double[] {ltvSchedule[iOlTV][ltv_list.size() - abs(idx_u) + 1], ltvSchedule[iOlTV][ltv_list.size() - abs(idx_u)]};
                            double[] ltv_bounds_d = new double[] {ltvSchedule[iOlTV][ltv_list.size() - abs(idx_d) + 1], ltvSchedule[iOlTV][ltv_list.size() - abs(idx_d)]};
                            double[] expectation_bounds_u_e = new double[2];
                            double[] expectation_bounds_d_e = new double[2];
                            double[] expectation_bounds_u_f = new double[2];
                            double[] expectation_bounds_d_f = new double[2];
                            double expected_u, expected_d;
                            for (int ie = 0; ie < ne; ie++){
                                for (int iFp = 0; iFp < nFICO; iFp++){
                                    for (int iep = 0; iep < ne; iep++){
                                        expectation_bounds_u_e[0] = expectation_bounds_u_e[0] + P_E.get(ie, iep) * V[age + 1][ia][iq][iOlTV][ltv_list.size() - abs(idx_u) + 1][iFp][iep];
                                        expectation_bounds_u_e[1] = expectation_bounds_u_e[1] + P_E.get(ie, iep) * V[age + 1][ia][iq][iOlTV][ltv_list.size() - abs(idx_u)][iFp][iep];
                                        expectation_bounds_d_e[0] = expectation_bounds_d_e[0] + P_E.get(ie, iep) * V[age + 1][ia][iq][iOlTV][ltv_list.size() - abs(idx_d) + 1][iFp][iep];
                                        expectation_bounds_d_e[1] = expectation_bounds_d_e[1] + P_E.get(ie, iep) * V[age + 1][ia][iq][iOlTV][ltv_list.size() - abs(idx_d)][iFp][iep];
                                    }
                                    expectation_bounds_u_f[0] = expectation_bounds_u_f[0] + P_F.get(iF, iFp) * expectation_bounds_u_e[0];
                                    expectation_bounds_u_f[1] = expectation_bounds_u_f[1] + P_F.get(iF, iFp) * expectation_bounds_u_e[1];
                                    expectation_bounds_d_f[0] = expectation_bounds_d_f[0] + P_F.get(iF, iFp) * expectation_bounds_d_e[0];
                                    expectation_bounds_d_f[1] = expectation_bounds_d_f[1] + P_F.get(iF, iFp) * expectation_bounds_d_e[1];
                                }
                                polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds_u, expectation_bounds_u_f);
                                expected_u = polynomialSplineFunction.value(ltv_u);
                                polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds_d, expectation_bounds_d_f);
                                expected_d = polynomialSplineFunction.value(ltv_d);
                                Expectation[age + 1][ia][iq][iOlTV][iLTV][iF][ie] = pi_u * expected_u + pi_d * expected_d;
                            }
                        }
                    }
                }
            }
        }
        return Expectation;
    }

    public static void main(String[] args){
        Grids grids = new Grids();
        double[][] creditSurface = grids.creditSurface("src/LifeCycleProblem/Data/CreditSurface50.csv", 7, 23);
        for (int f = 0; f < 7; f++){
            System.out.println(Arrays.toString(creditSurface[f]));
        }
    }
}
