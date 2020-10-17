package Housing;

import com.imsl.math.Spline2DInterpolate;
import sgpp.*;

import java.util.ArrayList;
import java.util.Arrays;

public class Control {

    Spline2DInterpolate creditSurface(String filePathForCS, int nFICO, int nOLTV){
        CreditSurface creditSurfaceInterpolation = new CreditSurface();
        double[] ficos = creditSurfaceInterpolation.ficos();
        double[] originationLTVs = creditSurfaceInterpolation.originationLTVs(filePathForCS, nOLTV);
        double[][] surface = creditSurfaceInterpolation.creditSurface(filePathForCS, nFICO, nOLTV);
        Spline2DInterpolate interpolate = new Spline2DInterpolate(ficos, originationLTVs, surface, 2, 2);
        return interpolate;
    }

    public static void main(String[] args){

        // credit surface
        String creditSurfaceData = "src/Housing/Data/CreditSurface.csv";
        int numberOfFicos = 7, numberOfOltvs = 24;
        Spline2DInterpolate creditSurface = new Control().creditSurface(creditSurfaceData, numberOfFicos, numberOfOltvs);

        // origination ltv grid
        double[] originationLtvGrid = new Grids().originationLtvGrid(creditSurfaceData, numberOfOltvs);

        // asset grid
        int na = 15;
        double amin = 1e-3, amax = 1;
        double[] assetGrid = new Grids().linspace(na, amin, amax);

        // housing quantity
        int nq = 15;
        double qmin = 1e-3, qmax = 1;
        double[] quantityGrid = new Grids().linspace(nq, qmin, qmax);

        // Tauchen
        int ne = 9;
        double sigma_eps = 0.02058, lambda_eps = 0.99, m = 1.5;
        double[] productivityGrid = new Tauchen().productivityGrid(ne, sigma_eps, lambda_eps, m);

        sgpp.LoadJSGPPLib.loadJSGPPLib();

        Terminal terminal = new Terminal();

        int dimension = 6;
        Grid gridTerminal = Grid.createLinearBoundaryGrid(dimension);
        GridStorage gridStorage = gridTerminal.getStorage();

        int level = 2;
        gridTerminal.getGenerator().regular(level);

        DataVector alphaTerminal = new DataVector(gridStorage.getSize());
        alphaTerminal.setAll(0.0);

        int age = 29;
        ArrayList<double[]> outputMatrix = new ArrayList<>();
        for (int step = 0; step < 1; step++){
            for (int i = 0; i < gridStorage.getSize(); i++) {
                GridPoint gp = gridStorage.getPoint(i);
                double[] outputTerminal = terminal.householdTerminal(age, gp.getStandardCoordinate(0), gp.getStandardCoordinate(1),
                                                                     gp.getStandardCoordinate(2), gp.getStandardCoordinate(3),
                                                                     gp.getStandardCoordinate(4), gp.getStandardCoordinate(5),
                                                                     creditSurface, assetGrid, productivityGrid);
                outputMatrix.add(outputTerminal);
                alphaTerminal.set(i, outputTerminal[1]);
            }
            jsgpp.createOperationHierarchisation(gridTerminal).doHierarchisation(alphaTerminal);
            SurplusRefinementFunctor functor = new SurplusRefinementFunctor(alphaTerminal, 1);
            gridTerminal.getGenerator().refine(functor);
            System.out.println("refinement step " + (step + 1) + ", new grid size: " + alphaTerminal.getSize());
            alphaTerminal.resizeZero(gridStorage.getSize());

        }

        new Grids().dataWriter(outputMatrix, "src/Housing/Data/OutputData.csv");

        //----------------------------------//
        age = age - 1;
        Grid gridFuture = gridTerminal;
        DataVector alphaFuture = alphaTerminal;

        Interim interim = new Interim();
        double[] outputInterim = interim.householdInterim(age, 0.5, 0.14, 0.60, 0, 0.6, 0,
                creditSurface, assetGrid, productivityGrid, quantityGrid, originationLtvGrid, gridFuture, alphaFuture);
        System.out.println(Arrays.toString(outputInterim));





    }

}


