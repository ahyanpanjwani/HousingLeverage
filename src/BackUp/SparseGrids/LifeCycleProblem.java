package Sparse;

import LifeCycleProblem.CartesianSet;
import org.apache.commons.lang3.ArrayUtils;
import sgpp.*;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class LifeCycleProblem {


    double[] xgrid(int nx, double xmin, double xmax){

        double[] xgrid = new double[nx];
        xgrid[0] = xmin;
        double xstep = (xmax - xmin) / (nx - 1);

        for (int i = 1; i < xgrid.length; i++){
            xgrid[i] = xgrid[i - 1] + xstep;
        }

        return xgrid;
    }

    // let's just do the terminal period first
    double value(double e, double x, int[] xgrid){

        double r = 0.03;
        double beta = 0.97;

        double wage = exp(e);
        double futureAsset;
        double bequestU;
        double bequestD;
        double expected;
        double consumption;
        double utility;
        double VV = pow(-10, 5);
        double CC = 0;

        for (int ixp = 0; ixp < xgrid.length; ixp++){
            futureAsset = xgrid[ixp] / 100d;
            bequestU = (1 + 0.04) * futureAsset;
            bequestD = (1 + 0.02) * futureAsset;
            expected = 0.5 * log(bequestU) + 0.5 * log(bequestD);
            consumption = (1 + r) * x + wage -  futureAsset;
            utility = log(consumption) + beta * expected;

            if (consumption <= 0){utility = pow(-10, 5);}
            if (utility >= VV){
                VV = utility;
                CC = consumption;
            }
        }

        return VV;
    }

    public static void main(String[] args){

        long startTime = System.nanoTime();

        LifeCycleProblem lifeCycleProblem = new LifeCycleProblem();


        sgpp.LoadJSGPPLib.loadJSGPPLib();


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
            output[z][2] = lifeCycleProblem.value(x0, x1, xgrid);
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
                alpha.set(g, lifeCycleProblem.value(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1), xgrid));
            });


            jsgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);
            SurplusRefinementFunctor functor = new SurplusRefinementFunctor(alpha, alpha.getSize(), 0.01);
            grid.getGenerator().refine(functor);
            System.out.println("new grid size: " + alpha.getSize());
            alpha.resizeZero(gridStorage.getSize());


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

        try(FileWriter fileWriter = new FileWriter("src/Sparse/Life.csv")){
            for (int z = 0; z < Z; z++){
                fileWriter.append(String.valueOf(exp(output[z][0])));
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
