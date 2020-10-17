package BackUp;

import Jama.Matrix;
import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.lang3.ArrayUtils;

import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

/**
 * This iteration, v.F., introduces precomputation of integrals, i.e. instead of computing integrals for expectations
 * on the go every time, expectations for the whole state space gets computed after each backward induction step
 * and then the relevant values are used in the computation of the optimization routine. The hypercube with the
 * expectations is calculated in the LifeCycleProblem.Grids.java file and called into LifeCycleProblem.LifeCycleProblemF.java file. The expectations are
 * purely over the productivity process (AR(1) discretized a la Tauchen (1986)) so the precomputation is possible.
 * There is no interpolation over LTV which I will introduce in the next iteration and precomputation may not be possible
 * there; we'll see.
 *
 * Moreover, I used this iteration to create the demand curve in DataWork/LeverageCycle/Update_29June2020.html using
 * the Walrasian logic of guessing a price -> create the credit surface -> find the quantity of housing demanded.
 * I ran this model for prices: 60 to 10 at increments of 5 units.
 * The output is in DataWork/LeverageCycle/DemandCurveSynced.csv; synced referencing the Walrasian logic.
 *
 * Currently, the file says "fix this", referring to the output produced by the backward induction. These are usually
 * the next period variables that need to be recorded. I will fix these in the next iteration of the model where
 * it will make more sense since we will be creating the expectations over price/LTV' anyway.
 *
 * As of July 10th 2020, this model is final in that I do not plan to make any further changes to this version of the
 * model at all. I will build on it in the next iteration of the model, LifeCycleProblem.LifeCycleProblemG.java.
 */

public class LifeCycleProblemF {

    // the main function for finding the value function given state and ambient variables
    double[] value_T(int age, int ia, int ie, int iq, int iCS, int in, int T, int na,
                   double[] agrid, double[] egrid, Matrix P, double[] qgrid, double[] lgrid, double[] mgrid, double[] ngrid,
                     double alpha, double beta, double w, double r, double pH, double mu, int N){

        double mortgage_payment;
        double balance;
        double ltv = 0;
        double balance_p;
        double ltv_p = 0;
        double utility;
        double consumption;
        double bequest;
        double VV = pow(-10, 5);
        double CC = 0;
        double QQ = 0;
        double AA = 0;

        // problem for households in the terminal period of life (choose risk-free assets to maximize bequest)
        if (age == T - 1){
            mortgage_payment = mgrid[iCS] * lgrid[iCS] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCS], - N));
            // balance and ltv in the next period (p = prime) i.e. the period after death when their assets will be liquidated for the bequest
            balance = lgrid[iCS] * pH * qgrid[iq] * (pow(1 + mgrid[iCS], N) - pow(1 + mgrid[iCS], ngrid[in])) / (pow(1 + mgrid[iCS], N) - 1);
            ltv = balance / (pH * qgrid[iq]);
            balance_p = lgrid[iCS] * pH * qgrid[iq] * (pow(1 + mgrid[iCS], N) - pow(1 + mgrid[iCS], ngrid[in] + 1)) / (pow(1 + mgrid[iCS], N) - 1);
            ltv_p = balance_p / (pH * mu * qgrid[iq]);

            for (int iap = 0; iap < na; iap++){
                consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - agrid[iap] - mortgage_payment;
                bequest = agrid[iap] + (pH * mu * qgrid[iq]) * (1 - ltv_p);

                utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * log(bequest);

                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV){
                    VV = utility;
                    CC = consumption;
                    QQ = qgrid[iq];
                    AA = agrid[iap];
                }
            }
        }

        double[] output = new double[15];
        output[0] = age;
        output[1] = agrid[ia];
        output[2] = egrid[ie];
        output[3] = qgrid[iq];
        output[4] = mgrid[iCS];
        output[5] = lgrid[iCS];
        output[6] = ltv;
        output[7] = VV;
        output[8] = CC;
        output[9] = AA;
        output[10] = QQ;
        output[11] = mgrid[iCS];
        output[12] = lgrid[iCS];
        output[13] = ltv_p;
        output[14] = 0d;
        return output;
    }

    double[] value_t(int age, int ia, int ie, int iq, int iCS, int in,
                     int T, int na, int nq, int nCS,
                     double[] agrid, double[] egrid, double[] qgrid, double[] lgrid, double[] mgrid, double[] ngrid,
                     double[][][][][][] E, double alpha, double beta, double w, double r, double pH, double mu, int N){

        double VV = pow(-10, 5);
        double CC = 0;
        double QQ = 0;
        double AA = 0;
        double ltv = 0;
        double ltv_p = 0;
        double status = 0;

        // problems for households in ther interim periods of life (current, refi, buy/sell)
        if (age < T - 1 && age > 0){

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

            double utility;
            double consumption;

            double expected_C;
            double expected_R;
            double expected_S;
            double expected_D;

            double mortgage_payment = mgrid[iCS] * lgrid[iCS] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCS], - N));
            double current_ltv = lgrid[iCS] * (pow(1 + mgrid[iCS], N) - pow(1 + mgrid[iCS], ngrid[in])) / (pow(1 + mgrid[iCS], N) - 1);

            for (int iap = 0; iap < na; iap++){
                //--------------STAY CURRENT: START--------------------//
                expected_C = E[age + 1][iap][iq][iCS][in][ie];
                consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - mortgage_payment - agrid[iap];
                utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected_C;
                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV_C){
                    VV_C = utility;
                    CC_C = consumption;
                    AA_C = agrid[iap];
                }
                //-----------------STAY CURRENT: END--------------------//

                //-----------------REFINANCE / BUY-SELL / DEFAULT: START---------------------//
                for (int iCSp = 0; iCSp < nCS; iCSp++){
                    expected_R = E[age + 1][iap][iq][iCSp][1][ie];
                    mortgage_payment = mgrid[iCSp] * lgrid[iCSp] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCSp], - N));
                    consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) + (lgrid[iCSp] - current_ltv) * pH * qgrid[iq] - agrid[iap] - mortgage_payment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected_R;
                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV_R){
                        VV_R = utility;
                        CC_R = consumption;
                        AA_R = agrid[iap];
                    }
                    //--------------------REFINANCE: END-----------------------------//


                    for (int iqp = 0; iqp < nq; iqp++){
                        //--------------------BUY/SELL: START----------------------------//
                        expected_S = E[age + 1][iap][iqp][iCSp][1][ie];
                        mortgage_payment = mgrid[iCSp] * lgrid[iCSp] * pH * qgrid[iqp] / (1 - pow(1 + mgrid[iCSp], - N));
                        consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) + (1 - current_ltv) * pH * qgrid[iq] - (1 - lgrid[iCSp]) * pH * qgrid[iqp] - agrid[iap] - mortgage_payment;
                        utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iqp]) + beta * expected_S;
                        if (consumption <= 0){utility = pow(-10, 5);}
                        if (utility >= VV_S){
                            VV_S = utility;
                            CC_S = consumption;
                            AA_S = agrid[iap];
                            QQ_S = qgrid[iqp];
                        }
                        //--------------------BUY/SELL: END----------------------------//

                        //--------------------DEFAULT: START---------------------------//
                        expected_D = E[age + 1][iap][iqp][iCSp][0][ie];
                        consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - (1 - lgrid[iCSp]) * pH * qgrid[iqp] - agrid[iap] - mortgage_payment;
                        utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected_D;
                        if (consumption <= 0){utility = pow(-10, 5);}
                        if (utility >= VV_D){
                            VV_D = utility;
                            CC_D = consumption;
                            AA_D = agrid[iap];
                            QQ_D = qgrid[iqp];
                        }
                        //--------------------DEFAULT: END---------------------------//
                    }
                }
                //-----------------REFINANCE / BUY-SELL: END---------------------//
            }

            // now find the option that provides maximum value (and corresponding policy function like consumption and quantity of housing owned)
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

        }

        double[] output = new double[15];
        output[0] = age;
        output[1] = agrid[ia];
        output[2] = egrid[ie];
        output[3] = qgrid[iq];
        output[4] = mgrid[iCS];
        output[5] = lgrid[iCS];
        output[6] = ltv;        // fix this
        output[7] = VV;
        output[8] = CC;
        output[9] = AA;
        output[10] = QQ;
        output[11] = mgrid[iCS]; // fix this
        output[12] = lgrid[iCS]; // fix this
        output[13] = ltv_p;      // fix this
        output[14] = status;
        return output;
    }


    // the routine for solving a new born household's problem.
    double[] value_0(int age, int ia, int ie,
                     int T, int na, int ne, int nq, int nCS, int nn,
                     double[] agrid, double[] egrid, Matrix P, double[] qgrid, double[] lgrid, double[] mgrid, double[] ngrid,
                     double[][][][][][] E, double alpha, double beta, double w, double r, double pH, double mu, int N){

        double expected;
        double mortage_payment;
        double downpayment;
        double consumption;
        double utility;
        double VV = pow(-10, 5);
        double CC = 0;
        double QQ = 0;
        double AA = 0;
        double MM = 0;
        double LL = 0;
        double ltv_p = 0;

        for (int iap = 0; iap < na; iap++){
            for (int iq = 0; iq < nq; iq++){
                for (int iCS = 0; iCS < nCS; iCS++){

                    expected = E[age + 1][iap][iq][iCS][1][ie];
                    mortage_payment = mgrid[iCS] * lgrid[iCS] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCS], - N));
                    downpayment = (1 - lgrid[iCS]) * pH * qgrid[iq];
                    consumption = agrid[ia] + w * exp(egrid[ie]) - agrid[iap] - mortage_payment - downpayment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected;

                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV){
                        VV = utility;
                        CC = consumption;
                        QQ = qgrid[iq];
                        AA = agrid[iap];
                        MM = mgrid[iCS];
                        LL = lgrid[iCS];
                    }

                }
            }
        }

        double[] output = new double[15];
        output[0] = age;
        output[1] = agrid[ia];
        output[2] = egrid[ie];
        output[3] = 0;
        output[4] = 0;
        output[5] = 0;
        output[6] = 0;
        output[7] = VV;
        output[8] = CC;
        output[9] = AA;
        output[10] = QQ;
        output[11] = MM;
        output[12] = LL;
        output[13] = ltv_p;      // fix this
        output[14] = 5;
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
        int pH = 60;           // house price (index)
        double muH = 1;             // house price growth rate (gross)

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
        double[] ngrid = grids.ngrid(N + 1, 0, N);               // mortgage age grid
        double[][][][][][] V = grids.V(T, na, ne, nq, nCS, N + 1);     // value function for rest of life
        double[][] V_0 = grids.V_0(na, ne);                                // value function for newborns
        double[][] Q_0 = grids.V_0(na, ne);                                // policy function for quantity bought

        // Building the state space
        Integer[][] matrix1 = new Integer[5][];
        matrix1[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix1[1] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        matrix1[2] = ArrayUtils.toObject(IntStream.range(0, nq).toArray());
        matrix1[3] = ArrayUtils.toObject(IntStream.range(0, nCS).toArray());
        matrix1[4] = ArrayUtils.toObject(IntStream.range(0, N).toArray());

        CartesianSet<Integer> stateSpace1 = new CartesianSet<>(matrix1);

        LifeCycleProblemF lifeCycleProblemF = new LifeCycleProblemF();


        int Z = (int) stateSpace1.getCount();
        double[][][] capture = new double[T][Z][15];

        // start the backward induction
        int age_T = T - 1;
        System.out.println("age: " + age_T);
        // run cohort-level value function iteration in parallel
        IntStream.range(0, Z).parallel().forEach(z ->{
            List<Integer> node = stateSpace1.get(z);
            V[age_T][node.get(0)][node.get(1)][node.get(2)][node.get(3)][node.get(4)] =
                    lifeCycleProblemF.value_T(age_T, node.get(0), node.get(1), node.get(2), node.get(3), node.get(4),
                    T, na, agrid, egrid, P, qgrid, lgrid, mgrid, ngrid, alpha, beta, w, r, pH, muH, N)[7];
        });

        double[][][][][][] Expected;

        for (int age = T - 2; age > 0; age--){
            System.out.println("age: " + age);
            Expected = grids.Expected(T, na, nq, nCS, N, ne, V, age, P);
            double[][][][][][] finalExpected_t = Expected;
            int age_t = age;
            ProgressBar.wrap(
            IntStream.range(0, Z).parallel(), "Task: ").forEach(z ->{
                List<Integer> node = stateSpace1.get(z);
                V[age_t][node.get(0)][node.get(1)][node.get(2)][node.get(3)][node.get(4)] =
                        lifeCycleProblemF.value_t(age_t, node.get(0), node.get(1), node.get(2), node.get(3), node.get(4),
                        T, na, nq, nCS,
                        agrid, egrid, qgrid, lgrid, mgrid, ngrid, finalExpected_t, alpha, beta, w, r, pH, muH, N)[7];
            });
        }

        // create another state space for the new born problem
        Integer[][] matrix0 = new Integer[2][];
        matrix0[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix0[1] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());

        CartesianSet<Integer> stateSpace0 = new CartesianSet<>(matrix0);

        // solve the ownership problem for newborns
        int age_0 = 0;
        System.out.println("age: " + age_0);
        double[][][][][][] finalExpected_0 = grids.Expected(T, na, nq, nCS, N, ne, V, age_0, P);
        IntStream.range(0, na * ne).parallel().forEach(z -> {
            List<Integer> node = stateSpace0.get(z);
            capture[age_0][z] = lifeCycleProblemF.value_0(age_0, node.get(0), node.get(1), T, na, ne, nq, nCS, N,
                    agrid, egrid, P, qgrid, lgrid, mgrid, ngrid, finalExpected_0, alpha, beta, w, r, pH, muH, N);
        });


        double Q_D = 0;
        for (int age = T - 1; age >= 0; age--){
            for (int z = 0; z < Z; z++){
                Q_D = Q_D + capture[age][z][10];
            }
        }
        System.out.println("total q demanded: " + Q_D);


        long endTime = System.nanoTime();
        System.out.println("run time = " + (endTime - startTime) * 1e-9 / 60 + " mins");
    }

}