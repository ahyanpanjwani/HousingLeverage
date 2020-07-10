import org.apache.commons.lang3.ArrayUtils;

import java.util.stream.IntStream;

public class Testing {

    public static void main(String[] args){

        int T = 10;     // lifespan (in years)
        int na = 15;    // number of asset grid points
        int ne = 9;     // number of productivity grid points
        int nq = 15;    // number of quantity grid points
        int nCS = 20;   // number of credit surface grid points
        int N = T - 1;  // mortgage term (same as lifespan but starting at zero, ending one period earlier)

        // Building the state space
        Integer[][] matrix1 = new Integer[5][];
        matrix1[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix1[1] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        matrix1[2] = ArrayUtils.toObject(IntStream.range(0, nq).toArray());
        matrix1[3] = ArrayUtils.toObject(IntStream.range(0, nCS).toArray());
        matrix1[4] = ArrayUtils.toObject(IntStream.range(0, N).toArray());

        CartesianSet<Integer> stateSpace1 = new CartesianSet<>(matrix1);
        int Z = (int) stateSpace1.getCount();
        for (int z = 0; z < Z; z++){
            System.out.println(z);
        }

    }
}
