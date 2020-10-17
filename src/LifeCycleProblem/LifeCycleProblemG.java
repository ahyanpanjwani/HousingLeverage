package LifeCycleProblem;

import Jama.Matrix;
import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.lang.Math.*;
import static java.lang.Math.max;

/**
 * In this version of the model, LifeCycleProblem.LifeCycleProblemG.java, I will introduce uncertainty over house prices in expectations.
 * For details on how to formally do this using a discretized Geometric Brownian Motion a la Binomial Tress, see
 * WriteUps/Expectations-and-Interpolation.nb.html and Modelling Derivatives in C++ (p.133-135).
 *
 * At this stage I will not use precomputation to enable this dimension of uncertainty; agents will calculate it on
 * the fly. Eventually I would like to have it done in precomputation before backward inducting at each time step.
 * We'll see.
 *
 * Not only do I wish to introduce uncertainty on price (and, hence, LTV) in this version, v.G., I would also like
 * to formally embed the piecewise linear interpolation mechanism for off-grid LTVs using binary search. Let's see how
 * this goes.
 *
 * This iteration of the model (v.G) does include price uncertainty and interpolation (piece-wise linear, binary search)
 * but the run time has increased significantly (10mins/cohort vs 40secs/cohort in v.F). After some testing, I realized
 * that the bottle neck is the interpolation step which I will try to move into the precomputation of integrals stage
 * instead of interpolating on the fly as it is here. I also had to allow for ltv > 1, this is now build into the
 * model.
 */


public class LifeCycleProblemG {

    // first the problem for the agents in their last period of life
    double[] value_T(int age, int ia, int ie, int iq, int iCS, int il, int T, int na,
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

    double[] value_t(int age, int ia, int ie, int iq, int iCS, int il,
                     int T, int na, int nq, int nCS,
                     double[] agrid, double[] egrid, double[] qgrid, double[] lgrid, double[] mgrid,
                     double[][][][][][] E, double[][] ltv_schedule, double alpha, double beta, double w, double r, double p_h, double mu, double sigma,
                     int N){

        // mortgage terms
        double mortgage_rate = mgrid[iCS];              // given the agent's point on the CS, back out the mortgage rate
        double origination_ltv = lgrid[iCS];            // given the agent's point on the CS, back out the origination ltv
        double quantity = qgrid[iq];                    // given the agent's point on the quantity grid, back out how much housing they hold
        double existing_mortgage_payment =              // given the current price of housing and how much housing the guy has, calculate the mortgage payment
                mortgage_rate * origination_ltv * p_h * quantity / (1 - pow(1 + mortgage_rate, - N));
        double current_ltv = ltv_schedule[iCS][il];     // given how many periods have elapsed ('in'), calculate the current ltv of the mortgage

        double ltv_u = 0, ltv_d = 0;                    // ltv in the up & down states i.e. the state of the world (in the next period) in which house prices go up or down
        double balance_p;                       // balance in the next period after making this period's payment, we will need this to calculate the ltv in the next period

        // in case of refi, sales, default
        double new_mortgage_rate, new_origination_ltv, new_mortgage_payment;
        // new_quantity is needed only in case of buy/sell and default
        double new_quantity;

        // price expectations terms (discretizing a GBM as a two-point process)
        double delta_x = sqrt(pow(sigma, 2) + pow(mu, 2));      // step for house price
        double pi_u = 0.5 + (0.5 * mu / delta_x);               // prob. of house price going up
        double pi_d = 0.5 - (0.5 * mu / delta_x);               // prob. of house price going down
        double p_h_u = p_h * exp(delta_x);                      // house price at the up node
        double p_h_d = p_h * exp(-delta_x);                     // house price at the down node

        // interpolation terms
        int idx_u, idx_d;
        double[] ltv_bounds;
        double[] expectation_bounds;
        LinearInterpolator linearInterpolator = new LinearInterpolator();
        PolynomialSplineFunction polynomialSplineFunction;

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

        double expected_C, expected_C_u, expected_C_d;
        double expected_R, expected_R_u, expected_R_d;
        double expected_S, expected_S_u, expected_S_d;
        double expected_D;


        /* this is the first main optimization loop: asset choice, other loops for choice variables will be
         * nested within this loop
         * */
        for (int iap = 0; iap < na; iap++){
            //#####################################################//
            //--------------STAY CURRENT: START--------------------//
            //#####################################################//

            /* first calculate the balance to be carried into the next period, having paid this period's mortgage
             * payment and the calculate the two possible ltvs tomorrow, depending on whether house prices go
             * up or down in the next period. The +/-1e-6 ensures that the expected ltv in the up and down
             * states are within the grid for the knife-edge cases. Ideally, the ltv's are off-the-grid
             * but we would have to extrapolate to find the expectation there but that can be risky.
             * Instead, I bound the ltv to a highest possible (and lowest) possible value. This makes sure that
             * the binary search function works properly and does not throw out-of-index exceptions. The
             * noise is small enough to not introduce massive errors (it kind of counters floating point errors
             * in fact)
             * */
            balance_p = (1 + mortgage_rate) * current_ltv * p_h * quantity - existing_mortgage_payment;
            ltv_u = min(ltv_schedule[iCS][0] - 1e-6, max(balance_p / (p_h_u * quantity), 0 + 1e-6));
            ltv_d = min(ltv_schedule[iCS][0] - 1e-6, max(balance_p / (p_h_d * quantity), 0 + 1e-6));

            //------------STAY CURRENT INTERPOLATION---------------//
            /* first, convert the array containing the ltv schedule to an array list
             * (this is needed to use the binary search library)
             * */
            List<Double> ltv_list = Arrays.stream(ltv_schedule[iCS]).boxed().collect(Collectors.toList());

            /* now sort the list in ascending order (originally, this ltv schedule
             * is in descending order; this is needed to apply the binary search function)
             * */
            Collections.sort(ltv_list);

            /* apply the binary search function. The function finds the ltv_upperbar > ltv_u.
             * ltv_lowerbar is simply the index for upperbar + 1
             * */
            idx_u = Collections.binarySearch(ltv_list, ltv_u);

            /* Having found the index, calculate the upper and lower bounds for ltv and corresponding
             * expectations
             * */
            ltv_bounds = new double[] {ltv_schedule[iCS][ltv_list.size() - abs(idx_u) + 1], ltv_schedule[iCS][ltv_list.size() - abs(idx_u)]};
            expectation_bounds = new double[] {E[age + 1][iap][iq][iCS][ltv_list.size() - abs(idx_u) + 1][ie], E[age + 1][iap][iq][iCS][ltv_list.size() - abs(idx_u)][ie]};

            /* Use the apache interpolation (linear) function to interpolate over the bounds
             * and then calculate the expectation based on the piece-wise linear interpolations
             * over the bounded interval
             * */
            polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds, expectation_bounds);
            expected_C_u = polynomialSplineFunction.value(ltv_u);


            /* now repeat exactly the same steps to calculate the expectations for the state in which
             * how prices will decline.
             * */
            idx_d = Collections.binarySearch(ltv_list, ltv_d);
            ltv_bounds = new double[] {ltv_schedule[iCS][ltv_list.size() - abs(idx_d) + 1], ltv_schedule[iCS][ltv_list.size() - abs(idx_d)]};
            expectation_bounds = new double[] {E[age + 1][iap][iq][iCS][ltv_list.size() - abs(idx_d) + 1][ie], E[age + 1][iap][iq][iCS][ltv_list.size() - abs(idx_d)][ie]};
            polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds, expectation_bounds);
            expected_C_d = polynomialSplineFunction.value(ltv_d);
            //------------END INTERPOLATION-------------------------//

            /* calculate the final expectation value as the weighted sum of the two expectations
             * in the two states (up and down) that I calculated using piece-wise linear interpolation
             * using binary search above.
             * */
            expected_C = pi_u * expected_C_u + pi_d * expected_C_d;

            // now calculate the consumption by re-arranging the budget constraint
            consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - existing_mortgage_payment - agrid[iap];

            /* now calculate the flow utility plus the expected continuation value if the agent sticks
             * with the existing mortgage.
             * */
            utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected_C;

            /* make sure that consumption is non-negative, if it is then assign a large negative
             * value so that it is off-equilibrium path.
             * */
            if (consumption <= 0){utility = pow(-10, 5);}

            /* see if this choice of assets is more valuable than the previous choice of assets
             * remember: this is all with in the asset-optimization loop
             * if the choice is better, then capture the value and important policy choices.
             * */
            if (utility >= VV_C){
                VV_C = utility;
                CC_C = consumption;
                AA_C = agrid[iap];
            }
            //######################################################//
            //-----------------STAY CURRENT: END--------------------//
            //######################################################//

            /* Start the nested loop for mortgage choice which will be used for refi, purchase
             * and default
             * */
            for (int iCSp = 0; iCSp < nCS; iCSp++){
                //#################################################//
                //----------------REFINANCE: START-----------------//
                //#################################################//

                new_mortgage_rate = mgrid[iCSp];
                new_origination_ltv = lgrid[iCSp];
                new_mortgage_payment = new_mortgage_rate * new_origination_ltv * p_h * quantity
                        / (1 - pow(1 + new_mortgage_rate, - N));
                balance_p = (1 + new_mortgage_rate) * new_origination_ltv * p_h * quantity - new_mortgage_payment;
                ltv_u = min(ltv_schedule[iCSp][0] - 1e-6, max(balance_p / (p_h_u * quantity), 0 + 1e-6));
                ltv_d = min(ltv_schedule[iCSp][0] - 1e-6, max(balance_p / (p_h_d * quantity), 0 + 1e-6));

                //-------------START: REFINANCE INTERPOLATION--------------//
                ltv_list = Arrays.stream(ltv_schedule[iCSp]).boxed().collect(Collectors.toList());
                Collections.sort(ltv_list);

                idx_u = Collections.binarySearch(ltv_list, ltv_u);
                ltv_bounds = new double[] {ltv_schedule[iCSp][ltv_list.size() - abs(idx_u) + 1], ltv_schedule[iCSp][ltv_list.size() - abs(idx_u)]};
                expectation_bounds = new double[] {E[age + 1][iap][iq][iCSp][ltv_list.size() - abs(idx_u) + 1][ie], E[age + 1][iap][iq][iCSp][ltv_list.size() - abs(idx_u)][ie]};
                polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds, expectation_bounds);
                expected_R_u = polynomialSplineFunction.value(ltv_u);

                idx_d = Collections.binarySearch(ltv_list, ltv_d);
                ltv_bounds = new double[] {ltv_schedule[iCSp][ltv_list.size() - abs(idx_d) + 1], ltv_schedule[iCSp][ltv_list.size() - abs(idx_d)]};
                expectation_bounds = new double[] {E[age + 1][iap][iq][iCSp][ltv_list.size() - abs(idx_d) + 1][ie], E[age + 1][iap][iq][iCSp][ltv_list.size() - abs(idx_d)][ie]};
                polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds, expectation_bounds);
                expected_R_d = polynomialSplineFunction.value(ltv_d);
                //-------------END: REFINANCE INTERPOLATION--------------//

                expected_R = pi_u * expected_R_u + pi_d * expected_R_d;
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
                    balance_p = (1 + new_mortgage_rate) * new_origination_ltv * p_h * new_quantity - new_mortgage_payment;
                    ltv_u = min(ltv_schedule[iCSp][0] - 1e-6, max(balance_p / (p_h_u * quantity), 0 + 1e-6));
                    ltv_d = min(ltv_schedule[iCSp][0] - 1e-6, max(balance_p / (p_h_d * quantity), 0 + 1e-6));

                    //-------------START: BUY/SELL INTERPOLATION--------------//
                    ltv_list = Arrays.stream(ltv_schedule[iCSp]).boxed().collect(Collectors.toList());
                    Collections.sort(ltv_list);

                    idx_u = Collections.binarySearch(ltv_list, ltv_u);
                    ltv_bounds = new double[] {ltv_schedule[iCSp][ltv_list.size() - abs(idx_u) + 1], ltv_schedule[iCSp][ltv_list.size() - abs(idx_u)]};
                    expectation_bounds = new double[] {E[age + 1][iap][iqp][iCSp][ltv_list.size() - abs(idx_u) + 1][ie], E[age + 1][iap][iqp][iCSp][ltv_list.size() - abs(idx_u)][ie]};
                    polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds, expectation_bounds);
                    expected_S_u = polynomialSplineFunction.value(ltv_u);

                    idx_d = Collections.binarySearch(ltv_list, ltv_d);
                    ltv_bounds = new double[] {ltv_schedule[iCSp][ltv_list.size() - abs(idx_d) + 1], ltv_schedule[iCSp][ltv_list.size() - abs(idx_d)]};
                    expectation_bounds = new double[] {E[age + 1][iap][iqp][iCSp][ltv_list.size() - abs(idx_d) + 1][ie], E[age + 1][iap][iqp][iCSp][ltv_list.size() - abs(idx_d)][ie]};
                    polynomialSplineFunction = linearInterpolator.interpolate(ltv_bounds, expectation_bounds);
                    expected_S_d = polynomialSplineFunction.value(ltv_d);
                    //-------------END: BUY/SELL INTERPOLATION--------------//

                    expected_S = pi_u * expected_S_u + pi_d * expected_S_d;
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
                /* the quantity of housing loop which contains the problem for buy/sell and
                 * default has ended (with the preceding brace).
                 * */
            }
            /* the credit surface loop which contains the refinance, buy/sell, and default problems
             * has ended (with the preceding brace)
             * */
        }
        /* the asset loop which contains all the four problems: stay current, refinance, buy/sell,
         * and default, has ended (with the preceding brace) i.e. all option values have been calculated
         * at this point. All that remains is to figure out the optimal choice and capture the numerical
         * values of all relevant/important variables (e.g. the value function and policy functions)
         * */

        /* now find the option that provides maximum value (and corresponding policy function
         * like consumption and quantity of housing owned). I simply use the built-in max function
         * to do so (the max function can only take two arguments so I nest it).
         * */

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

        int T = 10;     // lifespan (in years)
        int na = 15;    // number of asset grid points
        int ne = 9;     // number of productivity grid points
        int nq = 15;    // number of quantity grid points
        int nCS = 23;   // number of credit surface grid points
        int N = T - 1;  // mortgage term (same as lifespan but starting at zero, ending one period earlier)
        int nl = N + 5;    // number of contemporaneous LTV points (including some over shoot (the +5); LTV > 1)

        double amin = 0.1;      // minimum asset holding
        double amax = 4.0;      // maximum asset holding

        // Tauchen params
        double sigma_eps = 0.02058;
        double lambda_eps = 0.99;
        double m = 1.5;

        double qmin = 0.0001;      // minimum housing quantity
        double qmax = 0.3;       // maximum hosuing quantity

        // global params
        double alpha = 0.33;        // fraction of income spent on housing
        double beta = 0.97;         // discount factor
        double w = 1;               // wages
        double r = 0.03;            // risk free interest rate
        int pH = 60;                // house price (index)
        double muH = 0.00;             // house price drift
        double sigmaH = 0.06231;    // house price vol

        // file path for the credit surface
        String CSfilePath = "DataWork/LeverageCycle/StochCSRealistic" + String.valueOf(pH) + ".csv";


        // Calling all the grid routines
        Grids grids = new Grids();
        double[] agrid = grids.agrid(na, amin, amax);                     // asset grid
        double[] egrid = grids.egrid(ne, sigma_eps, lambda_eps, m);       // productivity grid
        Matrix P = grids.P(ne, sigma_eps, lambda_eps, m, egrid);          // productivity transition matrix
        double[] qgrid = grids.qgrid(nq, qmin, qmax);                     // quantity grid
        double[] lgrid = grids.lgrid(CSfilePath, nCS);                    // oltv grid
        double[] mgrid = grids.mgrid(CSfilePath, nCS);                    // mortgage rate grid
        double[][][][][][] V = grids.V(T, na, ne, nq, nCS, nl);            // value function for rest of life
        double[][][][][][] Expected;                                       // precomputation of integrals/expectations
        double[][] ltv_schedule = grids.LTV_Schedule(mgrid, lgrid, N, nCS);// for each oltv, the possible values for contemporaneous ltv with some overshoot i.e. ltv > 1


        // Building the state space
        Integer[][] matrix1 = new Integer[5][];
        matrix1[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix1[1] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        matrix1[2] = ArrayUtils.toObject(IntStream.range(0, nq).toArray());
        matrix1[3] = ArrayUtils.toObject(IntStream.range(0, nCS).toArray());
        matrix1[4] = ArrayUtils.toObject(IntStream.range(0, nl).toArray());

        CartesianSet<Integer> stateSpace1 = new CartesianSet<>(matrix1);

        LifeCycleProblemG lifeCycleProblemG = new LifeCycleProblemG();

        int Z = (int) stateSpace1.getCount();

        // first the problem for agents in their terminal period of life
        int age_T = T - 1;
        System.out.println("age: " + age_T);
        IntStream.range(0, Z).parallel().forEach(z -> {
            List<Integer> node = stateSpace1.get(z);
            V[age_T][node.get(0)][node.get(1)][node.get(2)][node.get(3)][node.get(4)] =
                    lifeCycleProblemG.value_T(age_T, node.get(0), node.get(1), node.get(2), node.get(3), node.get(4),
                            T, na, agrid, egrid, qgrid, lgrid, mgrid, ltv_schedule,
                            alpha, beta, w, r, pH, muH, sigmaH, N)[1];
        });

        for (int age = T - 2; age > 0; age--){
            System.out.println("age: " + age);
            Expected = grids.Expected(T, na, nq, nCS, nl, ne, V, age, P);
            int age_t = age;
            double[][][][][][] finalExpected = Expected;
            ProgressBar.wrap(
                    IntStream.range(0, Z).parallel(), "Task: ").forEach(z -> {
                List<Integer> node = stateSpace1.get(z);
                V[age_t][node.get(0)][node.get(1)][node.get(2)][node.get(3)][node.get(4)] =
                        lifeCycleProblemG.value_t(age_t, node.get(0), node.get(1), node.get(2), node.get(3), node.get(4),
                                T, na, nq, nCS,
                                agrid, egrid, qgrid, lgrid, mgrid, finalExpected, ltv_schedule, alpha, beta, w, r, pH, muH, sigmaH, N)[1];
            });
        }

        long endTime = System.nanoTime();
        System.out.println("run time = " + (endTime - startTime) * 1e-9 / 60 + " mins");

    }


}
