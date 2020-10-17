package BackUp.SparseGrids;

import BackUp.CartesianSet;
import org.apache.commons.lang3.ArrayUtils;
import sgpp.*;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.exp;

public class Optimization {

    double Household(double currentAsset, double productivity){

        sgpp.LoadJSGPPLib.loadJSGPPLib();
        sgpp.jsgpp.omp_set_num_threads(1);

        Objective objective = new Objective();
        objective.setCurrentAsset(currentAsset);
        objective.setProductivity(productivity);
        objective.setInterestRate(0.03);
        objective.setDiscountRate(0.97);

        OptNelderMead nelderMead = new OptNelderMead(objective, 100);
        nelderMead.optimize();
        //DataVector xOptNM = nelderMead.getOptimalPoint();
        final double fXOptNM = nelderMead.getOptimalValue();
        return -fXOptNM;
    }

    public static void main(String[] args){

        long startTime = System.nanoTime();

        Optimization optimization = new Optimization();

        sgpp.LoadJSGPPLib.loadJSGPPLib();

        int dim = 2;
        Grid grid = Grid.createLinearBoundaryGrid(dim);
        GridStorage gridStorage = grid.getStorage();

        int level = 5;
        grid.getGenerator().regular(level);
        System.out.println("Number of points (CSG): " + gridStorage.getSize());

        DataVector alpha = new DataVector(gridStorage.getSize());
        alpha.setAll(0.0);

        int G = (int) gridStorage.getSize();
        IntStream.range(0, G).parallel().forEach(g -> {
            GridPoint gp = gridStorage.getPoint(g);
            alpha.set(g, optimization.Household(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)));
        });

        jsgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);
        DataVector p = new DataVector(dim);
        OperationEval opEval = jsgpp.createOperationEval(grid);

        int N = 100;
        int[] xgrid = IntStream.range(0, N).toArray();
        Integer[][] matrix = new Integer[2][];
        matrix[0] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix[1] = ArrayUtils.toObject(IntStream.range(0, N).toArray());

        CartesianSet<Integer> stateSpace = new CartesianSet<>(matrix);
        int Z = (int) stateSpace.getCount();

        double[][] output = new double[Z][3];

        for (int z = 0; z < Z; z++){
            List<Integer> node = stateSpace.get(z);
            double x0 = xgrid[node.get(0)] / 100d;
            double x1 = xgrid[node.get(1)] / 100d;
            p.set(0, x0);
            p.set(1, x1);
            output[z][0] = x0;
            output[z][1] = exp(x1);
            output[z][2] = opEval.eval(alpha, p);
        }

        try(FileWriter fileWriter = new FileWriter("src/Sparse/LifeB.csv")){
            for (int z = 0; z < Z; z++){
                fileWriter.append(String.valueOf(output[z][0]));
                fileWriter.append(",");
                fileWriter.append(String.valueOf(output[z][1]));
                fileWriter.append(",");
                fileWriter.append(String.valueOf(output[z][2]));
                fileWriter.append("\n");
            }
        }catch (IOException ioe){ioe.printStackTrace();}

        long endTime = System.nanoTime();
        System.out.println("run time = " + (endTime - startTime) * 1e-9 / 60 + " mins");




    }
}
