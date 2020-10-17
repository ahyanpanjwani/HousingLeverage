package LifeCycleProblem;

/**
 * This iteration of the model (v.H) includes price uncertainty (via interpolation) through the ExpectedB function
 * in Grids.
 */

import Jama.Matrix;
import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.lang3.ArrayUtils;

import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class LifeCycleProblemH {

    double[] value_T(int age, int ia, int iq, int iCS, int il, int ie, int T, int na,
                     double[] agrid, double[] egrid, double[] qgrid, double[] lgrid, double[] mgrid,
                     double[][] ltv_schedule, double alpha, double beta, double w, double r,
                     double p_h, double mu, double sigma, int N){

        // mortgage terms
        double mortgage_rate;
        double origination_ltv;
        double quantity = 0;
        double mortgage_payment;
        double current_ltv;
        double balance_p;       // balance in the next period after making this period's payment
        double ltv_u, ltv_d;

        // price expectations terms
        double delta_x = sqrt(pow(sigma, 2) + pow(mu, 2));      // step for house price
        double pi_u = 0.5 + (0.5 * mu / delta_x);               // prob. of house price going up
        double pi_d = 0.5 - (0.5 * mu / delta_x);               // prob. of house price going down
        double p_h_u = p_h * exp(delta_x);                      // house price at the up node
        double p_h_d = p_h * exp(-delta_x);                     // house price at the down node

        // value function terms
        double utility;
        double consumption;
        double bequest_u;
        double bequest_d;
        double expected;
        double VV = pow(-10, 5);
        double CC = 0;
        double QQ = 0;
        double AA = 0;

        // problem for households in the terminal period of life (choose risk-free assets to maximize bequest)
        if (age == T - 1){

            // setting mortgage terms
            mortgage_rate = mgrid[iCS];
            origination_ltv = lgrid[iCS];
            quantity = qgrid[iq];
            mortgage_payment = mortgage_rate * origination_ltv * p_h * quantity / (1 - pow(1 + mortgage_rate, - N));
            current_ltv = ltv_schedule[iCS][il];

            // if the agent own the house outright, then there is no mortgage payment to make
            if (current_ltv == 0){mortgage_payment = 0;}

            /* balance and ltv in the next period (and states) after making a periodic payment in the current period
             * (max ensures this is not negative; needed for guys who will pay off mortgage in current period)
             */
            balance_p = (1 + mortgage_rate) * current_ltv * p_h * quantity - mortgage_payment;
            ltv_u = max(balance_p / (p_h_u * quantity), 0);
            ltv_d = max(balance_p / (p_h_d * quantity), 0);

            // start the optimization process over a'
            for (int iap = 0; iap < na; iap++){

                // re-writing the budget constraint
                consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - agrid[iap] - mortgage_payment;

                // bequest if house prices go up or down, respectively
                bequest_u = (1 + r) * agrid[iap] + (1 - ltv_u) * p_h_u * quantity;
                bequest_d = (1 + r) * agrid[iap] + (1 - ltv_d) * p_h_d * quantity;

                /* expectation of the log of bequest (this is important, we are taking expectation over the value of the
                 * bequest, not the value of the expected bequest i.e. E(log(b)) not log(E(b))
                 * (re. Jensen's inequality issues)
                 */
                expected = pi_u * log(bequest_u) + pi_d * log(bequest_d);

                // calculate the utility: flow utility (log/Cobb-Douglas) plus the expected log of bequest
                utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected;

                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV){
                    VV = utility;
                    CC = consumption;
                    QQ = qgrid[iq];
                    AA = agrid[iap];
                }
            }
        }

        double[] output = new double[5];
        output[0] = quantity;
        output[1] = VV;
        output[2] = CC;
        output[3] = AA;
        output[4] = QQ;
        return output;
    }

    double[] value_t(int age, int ia, int iq, int iCS, int il, int ie,
                     int na, int nq, int nCS,
                     double[] agrid, double[] qgrid, double[] mgrid, double[] lgrid, double[][] ltv_schedule,
                     double[] egrid, double[][][][][][] Expected,
                     double alpha, double beta, double w, double r, double p_h, int N){

        // mortgage terms
        double mortgage_rate = mgrid[iCS];              // given the agent's point on the CS, back out the mortgage rate
        double origination_ltv = lgrid[iCS];            // given the agent's point on the CS, back out the origination ltv
        double quantity = qgrid[iq];                    // given the agent's point on the quantity grid, back out how much housing they hold
        double existing_mortgage_payment =              // given the current price of housing and how much housing the guy has, calculate the mortgage payment
                mortgage_rate * origination_ltv * p_h * quantity / (1 - pow(1 + mortgage_rate, - N));
        double current_ltv = ltv_schedule[iCS][il];     // given how many periods have elapsed ('in'), calculate the current ltv of the mortgage

        // in case of refi, sales, default
        double new_mortgage_rate, new_origination_ltv, new_mortgage_payment;
        // new_quantity is needed only in case of buy/sell and default
        double new_quantity;

        // value function terms
        double utility;
        double consumption;
        double expected = 0;
        double default_penalty = 1;
        double VV = pow(-10, 5);

        // policy function terms
        double CC = 0;
        double QQ = 0;
        double AA = 0;
        double status = 0;

        // problems for households in their interim periods of life (current, refi, buy/sell)
        double VV_C = pow(-10, 5); //value from staying *C*urrent
        double CC_C = 0;
        double AA_C = 0;
        double VV_R = pow(-10, 5); //value from *R*efinancing
        double CC_R = 0;
        double AA_R = 0;
        double VV_S = pow(-10, 5); // value from *S*elling
        double CC_S = 0;
        double AA_S = 0;
        double QQ_S = 0;
        double VV_D = pow(-10, 5); // value from *D*efault
        double CC_D = 0;
        double AA_D = 0;
        double QQ_D = 0;

        double expected_C, expected_R, expected_S, expected_D;

        /* this is the first main optimization loop: asset choice, other loops for choice variables will be
         * nested within this loop
         * */
        for (int iap = 0; iap < na; iap++) {
            //#####################################################//
            //--------------STAY CURRENT: START--------------------//
            //#####################################################//
            expected_C = Expected[age + 1][iap][iq][iCS][il][ie];
            consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - existing_mortgage_payment - agrid[iap];
            utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected_C;
            if (consumption <= 0){utility = pow(-10, 5);}
            if (utility >= VV_C){
                VV_C = utility;
                CC_C = consumption;
                AA_C = agrid[iap];
            }

            for (int iCSp = 0; iCSp < nCS; iCSp++){
                //#################################################//
                //----------------REFINANCE: START-----------------//
                //#################################################//
                new_mortgage_rate = mgrid[iCSp];
                new_origination_ltv = lgrid[iCSp];
                new_mortgage_payment = new_mortgage_rate * new_origination_ltv * p_h * quantity
                        / (1 - pow(1 + new_mortgage_rate, - N));
                expected_R = Expected[age + 1][iap][iq][iCSp][4][ie];
                consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) + (new_origination_ltv - current_ltv) * p_h * quantity - agrid[iap] - new_mortgage_payment;
                utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected_R;
                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV_R){
                    VV_R = utility;
                    CC_R = consumption;
                    AA_R = agrid[iap];
                }
                //##################################################//
                //--------------------REFINANCE: END----------------//
                //##################################################//

                /* Start the nested loop for housing quantity. Buy/sell and and default will be
                 * nested within this loop
                 * */
                for (int iqp = 0; iqp < nq; iqp++){
                    //##################################################//
                    //---------------BUY/SELL: START--------------------//
                    //##################################################//
                    new_mortgage_rate = mgrid[iCSp];
                    new_origination_ltv = lgrid[iCSp];
                    new_quantity = qgrid[iqp];
                    new_mortgage_payment = new_mortgage_rate * new_origination_ltv * p_h * new_quantity
                            / (1 - pow(1 + new_mortgage_rate, - N));
                    expected_S = Expected[age + 1][iap][iq][iCSp][4][ie];
                    consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) + (1 - current_ltv) * p_h * quantity - (1 - new_origination_ltv) * p_h * new_quantity - agrid[iap] - new_mortgage_payment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iqp]) + beta * expected_S;
                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV_S){
                        VV_S = utility;
                        CC_S = consumption;
                        AA_S = agrid[iap];
                        QQ_S = qgrid[iqp];
                    }
                    //#####################################################//
                    //--------------------BUY/SELL: END--------------------//
                    //#####################################################//

                    //#####################################################//
                    //--------------------DEFAULT: START-------------------//
                    //#####################################################//

                    expected_D = expected_S;
                    consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - (1 - new_origination_ltv) * p_h * new_quantity - agrid[iap] - new_mortgage_payment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(quantity) - default_penalty + beta * expected_D;
                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV_D){
                        VV_D = utility;
                        CC_D = consumption;
                        AA_D = agrid[iap];
                        QQ_D = new_quantity;
                    }
                    //########################################################//
                    //--------------------DEFAULT: END------------------------//
                    //########################################################//
                }
            }
        }

        VV = max(max(VV_C, VV_S), max(VV_R, VV_D));
        if (VV == VV_C){
            CC = CC_C;
            AA = AA_C;
            QQ = qgrid[iq];
            status = 1;
        }else if (VV == VV_R){
            CC = CC_R;
            AA = AA_R;
            QQ = qgrid[iq];
            status = 2;
        }else if (VV == VV_S){
            CC = CC_S;
            AA = AA_S;
            QQ = QQ_S;
            status = 3;
        }else if (VV == VV_D){
            CC = CC_D;
            AA = AA_D;
            QQ = QQ_D;
            status = 4;
        }

        double[] output = new double[6];
        output[0] = quantity;
        output[1] = VV;
        output[2] = CC;
        output[3] = AA;
        output[4] = QQ;
        output[5] = status;
        return output;
    }

    public static void main(String[] args){
        long startTime = System.nanoTime();

        int T = 10;         // lifespan (in years)
        int na = 15;        // number of asset grid points
        int ne = 9;         // number of productivity grid points
        int nq = 15;        // number of quantity grid points
        int nCS = 23;       // number of credit surface grid points
        int N = T - 1;      // mortgage term (same as lifespan but starting at zero, ending one period earlier)
        int nl = N + 5;     // number of contemporaneous LTV points (including some over shoot (the +5); LTV > 1)

        double amin = 0.1;      // minimum asset holding
        double amax = 4.0;      // maximum asset holding

        // Tauchen params
        double sigma_eps = 0.02058;
        double lambda_eps = 0.99;
        double m = 1.5;

        double qmin = 0.0001;      // minimum housing quantity
        double qmax = 0.3;         // maximum hosuing quantity

        // global params
        double alpha = 0.33;        // fraction of income spent on housing
        double beta = 0.97;         // discount factor
        double w = 1;               // wages
        double r = 0.03;            // risk free interest rate
        int pH = 60;                // house price (index)
        double muH = 0.00;          // house price drift
        double sigmaH = 0.06231;    // house price vol

        // file path for the credit surface
        // use this file path for building JAR for cluster
        // String CSfilePath = "/home/ap2272/project/StochCSRealistic60.csv";
        // use this file path for running local
        String CSfilePath = "DataWork/LeverageCycle/StochCSRealistic" + String.valueOf(pH) + ".csv";

        Grids grids = new Grids();
        double[] agrid = grids.agrid(na, amin, amax);                      // asset grid
        double[] qgrid = grids.qgrid(nq, qmin, qmax);                      // quantity grid
        double[] mgrid = grids.mgrid(CSfilePath, nCS);                     // mortgage rate grid
        double[] lgrid = grids.lgrid(CSfilePath, nCS);                     // oltv grid
        double[][] ltv_schedule = grids.LTV_Schedule(mgrid, lgrid, N, nCS);// for each oltv, the possible values for contemporaneous ltv with some overshoot i.e. ltv > 1
        double[] egrid = grids.egrid(ne, sigma_eps, lambda_eps, m);        // productivity grid
        Matrix P = grids.P(ne, sigma_eps, lambda_eps, m, egrid);           // productivity transition matrix
        double[][][][][][] V = grids.V(T, na, nq, nCS, nl, ne);
        double[][][][][][] Expected;

        Integer[][] matrix1 = new Integer[5][];
        matrix1[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix1[1] = ArrayUtils.toObject(IntStream.range(0, nq).toArray());
        matrix1[2] = ArrayUtils.toObject(IntStream.range(0, nCS).toArray());
        matrix1[3] = ArrayUtils.toObject(IntStream.range(0, nl).toArray());
        matrix1[4] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());

        CartesianSet<Integer> stateSpace1 = new CartesianSet<>(matrix1);
        int Z = (int) stateSpace1.getCount();

        LifeCycleProblemH lifeCycleProblemH = new LifeCycleProblemH();

        // first the problem for agents in their terminal period of life
        int age_T = T - 1;
        System.out.println("age: " + age_T);
        IntStream.range(0, Z).parallel().forEach(z -> {
            List<Integer> node = stateSpace1.get(z);
            V[age_T][node.get(0)][node.get(1)][node.get(2)][node.get(3)][node.get(4)] =
                    lifeCycleProblemH.value_T(age_T, node.get(0), node.get(1), node.get(2), node.get(3), node.get(4),
                            T, na, agrid, egrid, qgrid, lgrid, mgrid, ltv_schedule,
                            alpha, beta, w, r, pH, muH, sigmaH, N)[1];
        });

        for (int age = T -2; age > 0; age--){
            //System.out.println("age: " + age);
            Expected = grids.ExpectedB(T, na, nq, nCS, nl, N, ne, V, age, P, pH, muH, sigmaH,
                    mgrid, lgrid, ltv_schedule);
            int age_t = age;
            double[][][][][][] finalExpected = Expected;
            ProgressBar.wrap(
                    IntStream.range(0, Z).parallel(), "age: " + age_t ).forEach(z -> {
                        List<Integer> node = stateSpace1.get(z);
                        V[age_t][node.get(0)][node.get(1)][node.get(2)][node.get(3)][node.get(4)] =
                        lifeCycleProblemH.value_t(age_t, node.get(0), node.get(1), node.get(2), node.get(3),
                                node.get(4), na, nq, nCS, agrid, qgrid, mgrid, lgrid,
                                ltv_schedule, egrid, finalExpected, alpha, beta, w, r, pH, N)[1];
            });
        }

        long endTime = System.nanoTime();
        System.out.println("run time = " + (endTime - startTime) * 1e-9 / 60 + " mins");

    }
}
