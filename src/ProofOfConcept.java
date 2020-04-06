import Jama.Matrix;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import java.lang.management.ManagementFactory;
import java.text.DecimalFormat;
import java.util.ArrayList;
import com.sun.management.OperatingSystemMXBean;


public class ProofOfConcept {

    public double ExcessDemand(double asset_endowment, double income, int fico_index, double price_guess,
                               double[] input_price_index, double[] input_yield_curve, double[] credit_surface){

        // SETTING UP THE CREDIT SURFACE

        double[] price_index = input_price_index;
        double[] yield_curve = input_yield_curve;

        double[] m_strip = credit_surface;
        ArrayList<Double> M = new ArrayList<>();
        for (int i = 0; i < m_strip.length; i++){
            M.add(i, m_strip[i]);
        }

        // SETTING UP THE REST OF THE STRUCTURE/GRIDS

        // Parameters
        double alpha = 0.33;
        double beta = 0.99;

        // Time Horizon
        int N = 30;

        // Prices
        double y = income; // fixed annual wage
        double r_f = 2; // risk free rate
        double p_0 = price_guess; // price guess
        double q_supply = 1; // fixed supply of housing, normalized to a unit

        // Grids
        double a_min = 1e-6;
        double a_max = 1.0;
        int na = 100;
        double[] asset = {a_min, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, a_max};
        ArrayList<Double> A = new ArrayList<>();
        for (int i = 0; i < asset.length; i++){
            A.add(i, asset[i]);
        }

        double a_0 = asset_endowment;  // original asset endowment at birth

        double LTV_min = 0.3;  // minimum possible origination LTV
        double LTV_max = 1.0;  // maximum possible origination LTV
        int nLTV = 36;         //  number of points on the LTV menu
        ArrayList<Double> LTV = new ArrayList<>();
        LTV.add(0, LTV_min);
        for (int i = 1; i < nLTV; i++){
            LTV.add(i, LTV.get(i - 1) + 0.02);
        }

        //double[] quality = {0.99, 0.995, 0.997, 0.998, 0.999, 1.0, 1.001, 1.002, 1.003, 1.005, 1.01};
        double[] quality = {0.0001, 0.00015, 0.0002, 0.00025, 0.0003, 0.00035, 0.0004, 0.00045, 0.0005, 0.00055, 0.0006,
                0.00065, 0.0007, 0.00075, 0.0008, 0.00085, 0.0009, 0.00095, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035};
        ArrayList<Double> Q = new ArrayList<>();
        for (int i = 0; i < quality.length; i++){
            Q.add(i, quality[i]);
        }

        // TERMINAL PERIOD t = N - 1 i.e. t = 29 years
        int t = N - 1;
        double[][][][] U_inter = new double[A.size()][Q.size()][M.size()][A.size()];
        double U_min = -1e6;

        for (int i = 0; i < A.size(); i++){
            for (int j = 0; j < Q.size(); j++){
                for (int k = 0; k < M.size(); k++){
                    for (int l = 0; l < A.size(); l++){
                        double a_current = A.get(i);
                        double q_0 = Q.get(j); // quantity of housing in hand
                        double m_0 = M.get(k) / 100; // mortgage rate upon origination
                        double ltv_0 = LTV.get(k); // LTV upon origination
                        double a_future = A.get(l);
                        double riskfree = r_f / 100; // annual terms
                        double mortgage_payment = (m_0 * ltv_0 * p_0 * q_0) / (1 - Math.pow((1 + m_0), -N));
                        double c = (1 + riskfree) * a_current + y + (p_0 * q_0) - a_future - mortgage_payment; // consumption in t=N
                        double b = a_future; // bequest
                        if (c > 0){
                            U_inter[i][j][k][l] = (1 - alpha) * Math.log(c) + alpha * Math.log(q_0) + beta * Math.log(b);
                        }else {
                            U_inter[i][j][k][l] = U_min;
                        }
                    }
                }
            }
        }

        double[][][] V_Terminal = new double[A.size()][Q.size()][M.size()];

        for (int i = 0; i < A.size(); i++){
            for (int j = 0; j < Q.size(); j++){
                for (int k = 0; k < M.size(); k++){
                    V_Terminal[i][j][k] = Doubles.max(U_inter[i][j][k]);
                }
            }
        }

        double[][][] V_future = new double[A.size()][Q.size()][M.size()];
        double[][][] V_current = new double[A.size()][Q.size()][M.size()];

        V_future = V_Terminal;

        // The main recursion for 0 < t < N
        while (t > 1){
            t--;
            //System.out.println("Period = " + t);
            for (int i = 0; i < A.size(); i++){
                for (int j = 0; j < Q.size(); j++){
                    for (int k = 0; k < M.size(); k++){
                        for (int l = 0; l < A.size(); l++){
                            double a_current = A.get(i);
                            double q_0 = Q.get(j);  // quantity of housing in hand
                            double m_0 = M.get(k) / 100;  // mortgage rate upon origination in annual terms
                            double ltv_0 = LTV.get(k);  // LTV upon origination
                            double a_future = A.get(l);
                            double riskfree = r_f / 100;  // annual terms
                            double mortgage_payment = (m_0 * ltv_0 * p_0 * q_0) / (1 - Math.pow((1 + m_0), -N)); //  # ANNUAL mortgage payment
                            double c = (1 + riskfree) * a_current + y - a_future - mortgage_payment;
                            if (c > 0){
                                U_inter[i][j][k][l] = (1 - alpha) * Math.log(c) + alpha * Math.log(q_0) + beta * V_future[l][j][k];
                            }else {
                                U_inter[i][j][k][l] = U_min;
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < A.size(); i++){
                for (int j = 0; j < Q.size(); j++){
                    for (int k = 0; k < M.size(); k++){
                        V_current[i][j][k] = Doubles.max(U_inter[i][j][k]);
                    }
                }
            }
            V_future = V_current;
        }

        // NOW THE ORIGINATION PERIOD, t = 0

        t = t - 1;
        //System.out.println("Period = " + t);

        double[][][] U_inter_origination = new double[Q.size()][M.size()][A.size()];
        for (int j = 0; j < Q.size(); j++){
            for (int k = 0; k < M.size(); k++){
                for (int l = 0; l < A.size(); l++){
                    double a_current = a_0; //fixed endowment at birth
                    double q_0 = Q.get(j); //qty of housing on hand
                    double m_0 = M.get(k) / 100; // mortgage rate upon origination in annual terms
                    double ltv_0 = LTV.get(k); // origination ltv
                    double a_future = A.get(l);
                    double riskfree = r_f / 100;
                    double mortgage_payment = (m_0 * ltv_0 * p_0 * q_0) / (1 - Math.pow((1 + m_0), -N)); // ANNUAL mortgage payment
                    double downpayment = (1 - ltv_0) * p_0 * q_0;
                    double c = a_current + y - a_future - mortgage_payment - downpayment;
                    if (c > 0){
                        U_inter_origination[j][k][l] = (1 - alpha) * Math.log(c) + alpha * Math.log(q_0) + beta * V_future[l][j][k];
                    }else {
                        U_inter_origination[j][k][l] = U_min;
                    }
                }
            }
        }

        //U_inter has the value for all (q, (m, LTV), a_1). I want to find the maximum possible value for a given (q, (m, LTV))

        double[][] U_q_m = new double[Q.size()][M.size()];

        for (int j = 0; j < Q.size(); j++){
            for (int k = 0; k < M.size(); k++){
                U_q_m[j][k] = Doubles.max(U_inter_origination[j][k]);
            }
        }

        // Now I want U(q) for a given q and the optimal mortgage, (m, LTV), for that q.
        double[] U_q = new double[Q.size()]; //store the welfare for each q \in Q
        int[] index_mortgage = new int[Q.size()]; //store the index of the optimal mortgage contract that maximizes U(q)

        // Now for every q, find the mortgage contract (m, LTV) that maximizes U(q)
        for (int j = 0; j < Q.size(); j++){
            U_q[j] = Doubles.max(U_q_m[j]);
            index_mortgage[j] = Doubles.indexOf(U_q_m[j], U_q[j]);
        }

        int optimal_q_index = Doubles.indexOf(U_q, Doubles.max(U_q));
        double q_demand = Q.get(optimal_q_index);

        //System.out.println("--------------");
        //System.out.println("a_0 = " + a_0 + ", y = " + y);
        //System.out.println("quantity demanded = " + q_demand);
        //System.out.println("excess demand = " + (q_demand - q_supply));
        //System.out.println("mkt clearing price = " + p_0);

        // Optimal mortgage contract
        double mortgage_rate = M.get(Ints.max(index_mortgage));
        double ltv = LTV.get(Ints.max(index_mortgage));

        DecimalFormat decimalFormat = new DecimalFormat("0.00");

        //System.out.println("mortgage rate = " + decimalFormat.format(mortgage_rate) + "%");
        //System.out.println("OLTV = " + decimalFormat.format(ltv));
        //System.out.println("--------------");

        double excess_demand = q_demand;

        return excess_demand;
    }

    public static void main(String args[]){

        long startTime = System.nanoTime();

        //double[] price_index = {108.792, 109.215, 109.643, 110.395, 111.248, 112.203, 113.274, 114.229, 114.991, 115.467, 115.682, 115.839, 116.056};
        double[] price_index = {0.515, 0.517, 0.519, 0.521, 0.523, 0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548};
        //double[] price_index = {0.517, 0.519, 0.521, 0.523, 0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548, 0.5595};
        //double[] price_index = {0.519, 0.521, 0.523, 0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548, 0.5595, 0.562};
        //double[] price_index = {0.521, 0.523, 0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548, 0.5595, 0.562, 0.5615};
        //double[] price_index = {0.523, 0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548, 0.5595, 0.562, 0.5615, 0.5605};
        //double[] price_index = {0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548, 0.5595, 0.562, 0.5615, 0.5605, 0.559};
        //double[] price_index = {0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548, 0.5595, 0.562, 0.5615, 0.5605, 0.559, 0.557};
        //double[] price_index = {0.529, 0.533, 0.536, 0.539, 0.544, 0.548, 0.5595, 0.562, 0.5615, 0.5605, 0.559, 0.557, 0.5555};
        //double[] price_index = {0.533, 0.536, 0.539, 0.544, 0.548, 0.5595, 0.562, 0.5615, 0.5605, 0.559, 0.557, 0.5555, 0.535};


        double[] yield_curve = {1.73, 1.74, 1.85, 2.28, 3.22, 3.75, 4.52, 4.97, 5.20, 5.86, 5.56};

        double endowment = 0.01;
        double income = 0.01;
        int fico_index = 6;

        String ficoFicoPath = "C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Data\\FicoFico.csv";

        Matrix creditSurface = new CreditSurfaceTestingB().CreditSurface("0", price_index, yield_curve, ficoFicoPath);
        double[][] cs = creditSurface.getArray();

        double excessDemand = new ProofOfConcept().ExcessDemand(endowment, income, fico_index, 0.104, price_index, yield_curve, cs[fico_index]);

        long endTime = System.nanoTime();

        long duration = (endTime - startTime);

        OperatingSystemMXBean bean = (com.sun.management.OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();

        System.out.println("CPU Usage = " + bean.getProcessCpuLoad() * 100 + "%");

        System.out.println(Math.round(duration * 1e-9)+ " seconds");
    }

}
