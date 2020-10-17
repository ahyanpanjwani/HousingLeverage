package Sparse;

import LifeCycleProblem.CartesianSet;
import org.apache.commons.lang3.ArrayUtils;
import sgpp.*;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class AdaptiveSG {

    private double f(double x0, double x1, double a, double b){
        return 1 / (abs(a - pow(x0, 4) - pow(x1, 4)) + b);
    }

    public static void main(String[] args){

        long startTime = System.nanoTime();

        AdaptiveSG adaptiveSG = new AdaptiveSG();

        double a = 0.5;
        double b = 0.1;

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
            output[z][2] = adaptiveSG.f(x0, x1, a, b);
        }

        System.out.println("Evaluating adaptive SG...");

        sgpp.LoadJSGPPLib.loadJSGPPLib();

        int dim = 2;
        Grid grid = Grid.createLinearBoundaryGrid(dim);
        GridStorage gridStorage = grid.getStorage();

        int level = 5;
        grid.getGenerator().regular(level);

        DataVector alpha = new DataVector(gridStorage.getSize());
        alpha.setAll(0.0);

        double epsilon = 0.01;
        double error = 1;

        while (error > epsilon){

            int G = (int) gridStorage.getSize();
            IntStream.range(0, G).parallel().forEach(g -> {
                GridPoint gp = gridStorage.getPoint(g);
                alpha.set(g, adaptiveSG.f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1), a, b));
            });


            jsgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);
            SurplusRefinementFunctor functor = new SurplusRefinementFunctor(alpha, alpha.getSize(), 0.01);
            grid.getGenerator().refine(functor);
            System.out.println("new grid size: " + alpha.getSize());
            alpha.resizeZero(gridStorage.getSize());

            OperationQuadrature operationQuadrature = jsgpp.createOperationQuadrature(grid);
            double res = operationQuadrature.doQuadrature(alpha);
            System.out.println("res_quadrature: " + res);

            DataVector p = new DataVector(dim);
            OperationEval opEval = jsgpp.createOperationEval(grid);

            for (int z = 0; z < Z; z++){
                List<Integer> node = stateSpace.get(z);
                p.set(0, xgrid[node.get(0)] / 100d);
                p.set(1, xgrid[node.get(1)] / 100d);
                output[z][3] = opEval.eval(alpha, p);
            }

            double sum = 0;
            for (int z = 0; z < Z; z++){
                sum = sum + pow(output[z][2] - output[z][3], 2);
            }
            error = sqrt(sum / Z);
            System.out.println("error: " + error);
        }

        try(FileWriter fileWriter = new FileWriter("src/Sparse/ASG.csv")){
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
