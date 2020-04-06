import org.javatuples.Triplet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Heterogeneity {

    public HashMap<Integer, Triplet<Double, Double, Integer>> Households(double[] asset_endowment_bins, double[] wages,
                                                                         int[] fico_index){

        int total_bins = asset_endowment_bins.length * wages.length * fico_index.length;

        List<Triplet<Double, Double, Integer>> H = new ArrayList<>();

        HashMap<Integer, Triplet<Double, Double, Integer>> Households = new HashMap<>();

        for (int i = 0; i < asset_endowment_bins.length; i++){
            for (int j = 0; j < wages.length; j++){
                for (int k = 0; k < fico_index.length; k++){
                    H.add(Triplet.with(asset_endowment_bins[i], wages[j], fico_index[k]));
                }
            }
        }

        for (int i = 0; i < total_bins; i++){
            Households.put(i, H.get(i));
        }

        return Households;
    }
}
