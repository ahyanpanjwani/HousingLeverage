import org.javatuples.Quartet;
import org.javatuples.Triplet;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.IntStream;

public class Control {

    public static void main(String args[]){

        long startTime = System.nanoTime();

        // Set up the grids/bins for heterogeneity

        double[] A = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
        double[] W = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
        int[] FICO_index = {0, 1, 2, 3, 4, 5, 6};

        HashMap<Integer, Triplet<Double, Double, Integer>> Households = new Heterogeneity().Households(A, W, FICO_index);

        double[] price_index = {0.515, 0.517, 0.519, 0.521, 0.523, 0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548};

        double[] yield_curve = {1.73, 1.74, 1.85, 2.28, 3.22, 3.75, 4.52, 4.97, 5.20, 5.86, 5.56};

        CreditSurface creditSurface = new CreditSurface();
        double[][] cs = creditSurface.CreditSurface("0", price_index, yield_curve).getArray();

        // This list will capture the (id, q_0, (m_0, LTV_0)) for all households
        List<Quartet<Integer, Double, Double, Double>> Demand = new ArrayList<>();

        BuyOnly buyOnly = new BuyOnly();

        double[] q_0 = new double[Households.size()];
        double[] m_0 = new double[Households.size()];
        double[] ltv_0 = new double[Households.size()];

        IntStream.range(0, Households.size()).parallel().forEach(z ->{
            Quartet<Integer, Double, Double, Double> demand = buyOnly.Demand(z, Households.get(z).getValue0(),
                    Households.get(z).getValue1(), cs[Households.get(z).getValue2()], 233.81,
                    price_index, yield_curve);
            q_0[z] = demand.getValue1();
            m_0[z] = demand.getValue2();
            ltv_0[z] = demand.getValue3();
        });

        DecimalFormat decimalFormat = new DecimalFormat("0.00");

        System.out.println("Total Demand = " + Arrays.stream(q_0).sum());
        System.out.println("Average mortgage rate = " + decimalFormat.format(Arrays.stream(m_0).average().getAsDouble()) + "%");
        System.out.println("Average OLTV = " + decimalFormat.format(Arrays.stream(ltv_0).average().getAsDouble()));

        long endTime = System.nanoTime();

        long duration = (endTime - startTime);

        System.out.println(Math.round(duration * 1e-9) + " seconds");
    }

}
