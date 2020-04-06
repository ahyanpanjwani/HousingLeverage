import Jama.Matrix;
import com.sun.management.OperatingSystemMXBean;
import org.javatuples.Triplet;

import java.lang.management.ManagementFactory;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.IntStream;

public class Control {

    public static void main(String args[]){

        long startTime = System.nanoTime();

        double[] A = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
        double[] W = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
        int[] FICO_index = {0, 1, 2, 3, 4, 5, 6};

        List<Triplet<Double, Double, Integer>> H = new ArrayList<Triplet<Double, Double, Integer>>();
        HashMap<Integer, Triplet<Double, Double, Integer>> Households = new HashMap<>();

        for (int i = 0; i < A.length; i++){
            for (int j = 0; j < W.length; j++){
                for (int k = 0; k < FICO_index.length; k++){
                    H.add(Triplet.with(A[i], W[j], FICO_index[k]));
                }
            }
        }

        for (int i = 0; i < A.length * W.length * FICO_index.length; i++){
            Households.put(i, H.get(i));
        }

        double[] price_index = {0.515, 0.517, 0.519, 0.521, 0.523, 0.525, 0.527, 0.529, 0.533, 0.536, 0.539, 0.544, 0.548};

        double[] yield_curve = {1.73, 1.74, 1.85, 2.28, 3.22, 3.75, 4.52, 4.97, 5.20, 5.86, 5.56};

        String ficoFicoPath = "C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Data\\FicoFico.csv";

        CreditSurface creditSurface = new CreditSurface();
        Matrix credit_surface = creditSurface.CreditSurface("0", price_index, yield_curve, ficoFicoPath);
        double[][] cs = credit_surface.getArray();

        double[] excess_demand = new double[A.length * W.length * FICO_index.length];

        IntStream.range(0, A.length * W.length * FICO_index.length).parallel().forEach(z ->{
            double excessDemand_h = new ProofOfConcept().ExcessDemand(Households.get(z).getValue0(),
                    Households.get(z).getValue1(), Households.get(z).getValue2(), 233.81, price_index, yield_curve,
                    cs[Households.get(z).getValue2()]);
            excess_demand[z] = excessDemand_h;
        });

        System.out.println("Max Individual Demand = " + Arrays.stream(excess_demand).max().getAsDouble());
        System.out.println("Min Individual Demand = " + Arrays.stream(excess_demand).min().getAsDouble());
        System.out.println("Total Demand = " + Arrays.stream(excess_demand).sum());

        OperatingSystemMXBean bean = (com.sun.management.OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();

        DecimalFormat decimalFormat = new DecimalFormat("0.00");

        System.out.println("CPU Usage = " + decimalFormat.format(bean.getProcessCpuLoad() * 100) + "%");

        long endTime = System.nanoTime();

        long duration = (endTime - startTime);

        System.out.println(Math.round(duration * 1e-9) + " seconds");
    }

}
