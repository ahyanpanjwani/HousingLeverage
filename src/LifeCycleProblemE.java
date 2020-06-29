import Jama.Matrix;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;
import static java.lang.Math.max;

/**
 * In this file, LifeCycleProblemE.java, I will update the D-version of the model to include default.
 * In v.D: stay current, refinance, buy/sell
 * In v.E: stay current, refinance, buy/sell, default
 * The option to default will be introduced as a nested loop with the buy/sell loop (under the loop for q).
 */

public class LifeCycleProblemE {

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

    // the main function for finding the value function given state and ambient variables
    double[] value(int age, int ia, int ie, int iq, int iCS, int in,
                   int T, int na, int ne, int nq, int nCS, int nn,
                   double[] agrid, double[] egrid, Matrix P, double[] qgrid, double[] lgrid, double[] mgrid, double[] ngrid,
                   double[][][][][][] V, double alpha, double beta, double w, double r, double pH, double mu, int N){

        double expected;
        double utility;
        double consumption;
        double bequest;
        double VV = pow(-10, 5);
        double CC = 0;
        double QQ = 0;
        double BB = 0;

        // problem for households in the terminal period of life (choose risk-free assets to maximize bequest)
        if (age == T - 1){
            double mortgage_payment = mgrid[iCS] * lgrid[iCS] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCS], - N));
            // balance and ltv in the next period (p = prime) i.e. the period after death when their assets will be liquidated for the bequest
            double balancep = lgrid[iCS] * pH * qgrid[iq] * (pow(1 + mgrid[iCS], N) - pow(1 + mgrid[iCS], ngrid[in])) / (pow(1 + mgrid[iCS], N) - 1);
            double ltvp = balancep / (pH * mu * qgrid[iq]);

            for (int iap = 0; iap < na; iap++){
                expected = 0;
                consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - agrid[iap] - mortgage_payment;
                bequest = agrid[iap] + (pH * mu * qgrid[iq]) * (1 - ltvp);

                utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * log(bequest);

                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV){
                    VV = utility;
                    CC = consumption;
                    BB = bequest;
                }
            }
        }

        double VV_C = pow(-10, 5); //value from staying *C*urrent
        double CC_C = 0;
        double VV_R = pow(-10, 5); //value from *R*efinancing
        double CC_R = 0;
        double VV_S = pow(-10, 5); // value from *S*elling
        double CC_S = 0;
        double QQ_S = 0;
        double VV_D = pow(-10, 5); // value from *D*efault
        double CC_D = 0;
        double QQ_D = 0;

        // problems for households in ther interim periods of life (current, refi, buy/sell)
        if (age < T - 1 && age >= 1){
            double mortgage_payment = mgrid[iCS] * lgrid[iCS] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCS], - N));
            double current_balance = lgrid[iCS] * pH * qgrid[iq] * (pow(1 + mgrid[iCS], N) - pow(1 + mgrid[iCS], ngrid[in])) / (pow(1 + mgrid[iCS], N) - 1);
            double current_ltv = current_balance / (pH * qgrid[iq]);

            for (int iap = 0; iap < na; iap++){
                //--------------STAY CURRENT: START--------------------//
                expected = 0;
                for (int iep = 0; iep < ne; iep++){
                    expected = expected + P.get(ie, iep) * V[age + 1][iap][iep][iq][iCS][in + 1];
                }

                consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - mortgage_payment - agrid[iap];
                utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected;
                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV_C){
                    VV_C = utility;
                    CC_C = consumption;
                }
                //-----------------STAY CURRENT: END--------------------//

                //-----------------REFINANCE / BUY-SELL / DEFAULT: START---------------------//
                for (int iCSp = 0; iCSp < nCS; iCSp++){
                    expected = 0;
                    //------------------REFINANCE: START---------------------------//
                    for (int iep = 0; iep < ne; iep++){
                        expected = expected + P.get(ie, iep) * V[age + 1][iap][iep][iq][iCSp][1];
                    }
                    mortgage_payment = mgrid[iCSp] * lgrid[iCSp] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCSp], - N));
                    consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) + (lgrid[iCSp] - current_ltv) * pH * qgrid[iq] - agrid[iap] - mortgage_payment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected;
                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV_R){
                        VV_R = utility;
                        CC_R = consumption;
                    }
                    //--------------------REFINANCE: END-----------------------------//


                    for (int iqp = 0; iqp < nq; iqp++){
                        expected = 0;
                        for (int iep = 0; iep < ne; iep++){
                            expected = expected + P.get(ie, iep) * V[age + 1][iap][iep][iqp][iCSp][1];
                        }
                        //--------------------BUY/SELL: START----------------------------//
                        mortgage_payment = mgrid[iCSp] * lgrid[iCSp] * pH * qgrid[iqp] / (1 - pow(1 + mgrid[iCSp], - N));
                        consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) + (1 - current_ltv) * pH * qgrid[iq] - (1 - lgrid[iCSp]) * pH * qgrid[iqp] - agrid[iap] - mortgage_payment;
                        utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iqp]) + beta * expected;
                        if (consumption <= 0){utility = pow(-10, 5);}
                        if (utility >= VV_S){
                            VV_S = utility;
                            CC_S = consumption;
                            QQ_S = qgrid[iqp];
                        }
                        //--------------------BUY/SELL: END----------------------------//

                        //--------------------DEFAULT: START---------------------------//
                        consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - (1 - lgrid[iCSp]) * pH * qgrid[iqp] - agrid[iap] - mortgage_payment;
                        utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected;
                        if (consumption <= 0){utility = pow(-10, 5);}
                        if (utility >= VV_D){
                            VV_D = utility;
                            CC_D = consumption;
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
                QQ = qgrid[iq];
                System.out.println("curr");
            }else if (VV == VV_R){
                CC = CC_R;
                QQ = qgrid[iq];
                System.out.println("refi");
            }else if (VV == VV_S){
                CC = CC_S;
                QQ = QQ_S;
                System.out.println("sell" + " , " + " q: " + qgrid[iq] +  " , q': " + QQ);
            }else if (VV == VV_D){
                CC = CC_D;
                QQ = QQ_D;
                System.out.println("deft");
            }

        }


        double[] output = new double[3];
        output[0] = VV;
        output[1] = CC;
        output[2] = QQ;
        return output;
    }

    // the routine for solving a new born household's problem.
    double[] value_0(int age, int ia, int ie,
                     int T, int na, int ne, int nq, int nCS, int nn,
                     double[] agrid, double[] egrid, Matrix P, double[] qgrid, double[] lgrid, double[] mgrid, double[] ngrid,
                     double[][][][][][] V, double alpha, double beta, double w, double r, double pH, double mu, int N){

        double expected;
        double mortage_payment;
        double downpayment;
        double consumption;
        double utility;
        double VV = pow(-10, 5);
        double CC = 0;
        double QQ = 0;

        for (int iap = 0; iap < na; iap++){
            for (int iq = 0; iq < nq; iq++){
                for (int iCS = 0; iCS < nCS; iCS++){
                    expected = 0;
                    for (int iep = 0; iep < ne; iep++){
                        expected = expected + P.get(ie, iep) * V[age + 1][iap][iep][iq][iCS][1];
                    }

                    mortage_payment = mgrid[iCS] * lgrid[iCS] * pH * qgrid[iq] / (1 - pow(1 + mgrid[iCS], - N));
                    downpayment = (1 - lgrid[iCS]) * pH * qgrid[iq];
                    consumption = agrid[ia] + w * exp(egrid[ie]) - agrid[iap] - mortage_payment - downpayment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(qgrid[iq]) + beta * expected;

                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV){
                        VV = utility;
                        CC = consumption;
                        QQ = qgrid[iq];
                    }

                }
            }
        }

        double[] out = new double[3];
        out[0] = VV;
        out[1] = CC;
        out[2] = QQ;
        return out;
    }

    public static void main(String args[]){

        long startTime = System.nanoTime();

        int T = 10;     // lifespan (in years)
        int na = 20;    // number of asset grid points
        int ne = 9;     // number of productivity grid points
        int nq = 20;    // number of quantity grid points
        int nCS = 26;   // number of credit surface grid points
        int N = T - 1;  // mortgage term (same as lifespan but starting at zero, ending one period earlier)

        double amin = 0.1;      // minimum asset holding
        double amax = 4.0;      // maximum asset holding

        // Tauchen params
        double sigma_eps = 0.02058;
        double lambda_eps = 0.99;
        double m = 1.5;

        double qmin = 0.0001;      // minimum housing quantity
        double qmax = 0.3;       // maximum hosuing quantity

        // file path for the credit surface
        String CSfilePath = "DataWork/LeverageCycle/StochCSRealistic.csv";

        // global params
        double alpha = 0.33;        // fraction of income spent on housing
        double beta = 0.97;         // discount factor
        double w = 1;               // wages
        double r = 0.03;            // risk free interest rate
        double pH = 30;           // house price (index)
        double muH = 1;             // house price growth rate (gross)


        // Calling all the grid routines
        LifeCycleProblemE lifeCycleProblemE= new LifeCycleProblemE();
        double[] agrid = lifeCycleProblemE.agrid(na, amin, amax);                     // asset grid
        double[] egrid = lifeCycleProblemE.egrid(ne, sigma_eps, lambda_eps, m);       // productivity grid
        Matrix P = lifeCycleProblemE.P(ne, sigma_eps, lambda_eps, m, egrid);          // productivity transition matrix
        double[] qgrid = lifeCycleProblemE.qgrid(nq, qmin, qmax);                     // quantity grid
        double[] lgrid = lifeCycleProblemE.lgrid(CSfilePath, nCS);                    // oltv grid
        double[] mgrid = lifeCycleProblemE.mgrid(CSfilePath, nCS);                    // mortgage rate grid
        double[] ngrid = lifeCycleProblemE.ngrid(N + 1, 0, N);               // mortgage age grid
        double[][][][][][] V = lifeCycleProblemE.V(T, na, ne, nq, nCS, N + 1);     // value function for rest of life
        double[][] V_0 = lifeCycleProblemE.V_0(na, ne);                                // value function for newborns
        double[][] Q_0 = lifeCycleProblemE.V_0(na, ne);                                // policy function for quantity bought

        // Building the stat space
        Integer[][] matrix1 = new Integer[5][];
        matrix1[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix1[1] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        matrix1[2] = ArrayUtils.toObject(IntStream.range(0, nq).toArray());
        matrix1[3] = ArrayUtils.toObject(IntStream.range(0, nCS).toArray());
        matrix1[4] = ArrayUtils.toObject(IntStream.range(0, N).toArray());

        CartesianSet<Integer> stateSpace1 = new CartesianSet<>(matrix1);


        // start the backward induction
        for (int age = T - 1; age >= 1; age--){
            int Z = (int) stateSpace1.getCount();
            int finalAge = age;
            System.out.println("age = " + age);
            // run cohort-level value function iteration in parallel
            IntStream.range(0, Z).parallel().forEach(z ->{
                List<Integer> node = stateSpace1.get(z);
                V[finalAge][node.get(0)][node.get(1)][node.get(2)][node.get(3)][node.get(4)]
                        = lifeCycleProblemE.value(finalAge, node.get(0), node.get(1), node.get(2), node.get(3), node.get(4),
                        T, na, ne, nq, nCS, N,
                        agrid, egrid, P, qgrid, lgrid, mgrid, ngrid, V, alpha, beta, w, r, pH, muH, N)[0];
            });

        }

        // create another state space for the new born problem
        Integer[][] matrix0 = new Integer[2][];
        matrix0[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix0[1] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());

        CartesianSet<Integer> stateSpace2 = new CartesianSet<>(matrix0);

        // solve the ownership problem for newborns
        int age = 0;
        System.out.println("age = " + age);
        IntStream.range(0, na * ne).parallel().forEach(z -> {
            List<Integer> node = stateSpace2.get(z);
            V_0[node.get(0)][node.get(1)] = lifeCycleProblemE.value_0(age, node.get(0), node.get(1), T, na, ne, nq, nCS, N,
                    agrid, egrid, P, qgrid, lgrid, mgrid, ngrid, V, alpha, beta, w, r, pH, muH, N)[0];
            Q_0[node.get(0)][node.get(1)] = lifeCycleProblemE.value_0(age, node.get(0), node.get(1), T, na, ne, nq, nCS, N,
                    agrid, egrid, P, qgrid, lgrid, mgrid, ngrid, V, alpha, beta, w, r, pH, muH, N)[2];
            //System.out.println(Q_0[node.get(0)][node.get(1)] + " , " + agrid[node.get(0)] + " , " + egrid[node.get(1)]);
        });

        long endTime = System.nanoTime();
        System.out.println("run time = " + (endTime - startTime) * 1e-9 / 60 + " mins");


    }


}
