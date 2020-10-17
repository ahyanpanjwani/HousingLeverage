package Sparse;

import LifeCycleProblem.CartesianSet;
import org.apache.commons.lang3.ArrayUtils;
import sgpp.*;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class ClassicalSG {

    private static double f(double x0, double x1){
        return 1 / (abs(0.5 - pow(x0, 4) - pow(x1, 4)) + 0.1);
        //return 16.0 * (x0 - 1.0) * x0 * (x1 - 1.0) * x1;
    }

    public static void main(String[] args){
        long startTime = System.nanoTime();

        int N = 100;
        int[] xgrid = IntStream.range(0, N).toArray();
        Integer[][] matrix = new Integer[2][];
        matrix[0] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix[1] = ArrayUtils.toObject(IntStream.range(0, N).toArray());

        CartesianSet<Integer> stateSpace = new CartesianSet<>(matrix);
        int Z = (int) stateSpace.getCount();

        double[][] output = new double[Z][4];

        System.out.println("Evaluating Cartesian...");
        double x0, x1;
        for (int z = 0; z < Z; z++){
            List<Integer> node = stateSpace.get(z);
            x0 = xgrid[node.get(0)] / 100d;
            x1 = xgrid[node.get(1)] / 100d;
            output[z][0] = x0;
            output[z][1] = x1;
            output[z][2] = f(x0, x1);
        }

        System.out.println("Evaluating classical SG...");

        sgpp.LoadJSGPPLib.loadJSGPPLib();

        int dim = 2;
        Grid grid = Grid.createLinearBoundaryGrid(dim);
        GridStorage gridStorage = grid.getStorage();

        int level = 14;
        grid.getGenerator().regular(level);
        System.out.println("Number of points (CSG): " + gridStorage.getSize());

        DataVector alpha = new DataVector(gridStorage.getSize());
        alpha.setAll(0.0);

        /*
        for (int i = 0; i < gridStorage.getSize(); i++) {
            GridPoint gp = gridStorage.getPoint(i);
            alpha.set(i, f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)));
        }
        */


        int G = (int) gridStorage.getSize();
        IntStream.range(0, G).parallel().forEach(g -> {
            GridPoint gp = gridStorage.getPoint(g);
            alpha.set(g, f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)));
        });


        jsgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);
        DataVector p = new DataVector(dim);
        OperationEval opEval = jsgpp.createOperationEval(grid);

        OperationQuadrature operationQuadrature = jsgpp.createOperationQuadrature(grid);
        double res = operationQuadrature.doQuadrature(alpha);
        System.out.println("res_quadrature: " + res);
        

        for (int z = 0; z < Z; z++){
            List<Integer> node = stateSpace.get(z);
            p.set(0, xgrid[node.get(0)] / 100d);
            p.set(1, xgrid[node.get(1)] / 100d);
            output[z][3] = opEval.eval(alpha, p);
        }


        double error = 0;
        for (int z = 0; z < Z; z++){
            error = error + pow(output[z][2] - output[z][3], 2);
        }
        error = sqrt(error / Z);
        System.out.println("error (CSG): " + error);

        try(FileWriter fileWriter = new FileWriter("src/Sparse/CSG.csv")){
            for (int z = 0; z < Z; z++){
                fileWriter.append(String.valueOf(output[z][0]));
                fileWriter.append(",");
                fileWriter.append(String.valueOf(output[z][1]));
                fileWriter.append(",");
                fileWriter.append(String.valueOf(output[z][2]));
                fileWriter.append(",");
                fileWriter.append(String.valueOf(output[z][3]));
                fileWriter.append("\n");
            }
        }catch (IOException ioe){ioe.printStackTrace();}

        long endTime = System.nanoTime();
        System.out.println("run time = " + (endTime - startTime) * 1e-9 / 60 + " mins");

    }
}
