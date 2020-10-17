package Sparse;

import Housing.CreditSurface;
import com.imsl.math.Spline2DInterpolate;
import sgpp.*;

import java.util.ArrayList;

public class Control {

    Spline2DInterpolate creditSurface(String filePathForCS, int nFICO, int nOLTV){
        Housing.CreditSurface creditSurfaceInterpolation = new CreditSurface();
        double[] ficos = creditSurfaceInterpolation.ficos();
        double[] originationLTVs = creditSurfaceInterpolation.originationLTVs(filePathForCS, nOLTV);
        double[][] surface = creditSurfaceInterpolation.creditSurface(filePathForCS, nFICO, nOLTV);
        Spline2DInterpolate interpolate = new Spline2DInterpolate(ficos, originationLTVs, surface, 2, 2);
        return interpolate;
    }

    double[] terminal(int age, double currentAsset, double currentQuantitiy, double currentOriginationLTV,
                  double ficoScore, double currentLTV, double currentProductivity,
                  Spline2DInterpolate creditSurface, double[] productivityGrid){

        ObjectiveFunctionTerminal householdTerminal = new ObjectiveFunctionTerminal();
        householdTerminal.setAge(age);
        householdTerminal.setCurrentAsset(currentAsset);
        householdTerminal.setCurrentQuantity(currentQuantitiy);
        householdTerminal.setCurrentOriginationLTV(currentOriginationLTV);
        householdTerminal.setFicoScore(ficoScore);
        householdTerminal.setCurrentLTV(currentLTV);
        householdTerminal.setCurrentProductivity(currentProductivity);
        householdTerminal.setCreditSurface(creditSurface);
        householdTerminal.setProductivityGrid(productivityGrid);

        householdTerminal.housePrice = 40;

        // dimension of the domain
        final long dim = householdTerminal.getNumberOfParameters();
        // B-spline degree
        final long p = 3;
        // max number of grid points
        final long N = 60;
        // adaptivity of grid generation
        final double gamma = 0.95;

        // First, we define a grid with modified B-spline basis functions and an iterative grid generator, which can generate the grid adaptively.
        Grid grid = Grid.createModBsplineGrid(dim, p);
        OptIterativeGridGeneratorRitterNovak gridGen = new OptIterativeGridGeneratorRitterNovak(householdTerminal, grid, N, gamma);

        // With the iterative grid generator, we generate adaptively a sparse grid.
        gridGen.generate();

        // Then, we hierarchize the function values to get hierarchical B-spline coefficients of the B-spline sparse grid interpolant \hat{f} : [0, 1]^d -> R
        final DataVector functionValues = gridGen.getFunctionValues();
        DataVector coeffs = new DataVector(functionValues.getSize());
        HierarchisationSLE hierarchisationSLE = new HierarchisationSLE(grid);
        AutoSLESolver sleSolver = new AutoSLESolver();

        sleSolver.solve(hierarchisationSLE, gridGen.getFunctionValues(), coeffs);

        // Now define the interpolant $\hat{f}$ and its gradient, $\nabla \hat{f}$ for use with the gradient method (steepest descent).
        InterpolantScalarFunction householdTerminalHat = new InterpolantScalarFunction(grid, coeffs);
        InterpolantScalarFunctionGradient householdTerminalHatGradient = new InterpolantScalarFunctionGradient(grid, coeffs);
        OptGradientDescent gradientDescent = new OptGradientDescent(householdTerminalHat, householdTerminalHatGradient);
        DataVector x0 = new DataVector(dim);
        double fX0;

        GridStorage gridStorage = grid.getStorage();
        // index of grid point with minimal function value
        int x0Index = 0;
        fX0 = functionValues.get(0);
        for (int i = 1; i < functionValues.getSize(); i++) {
            if (functionValues.get(i) < fX0) {
                fX0 = functionValues.get(i);
                x0Index = i;
            }
        }
        x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));

        // We apply the gradient method and print the results.
        gradientDescent.setStartingPoint(x0);
        gradientDescent.optimize();
        DataVector xOpt = gradientDescent.getOptimalPoint();
        final double value = - householdTerminal.eval(xOpt);
        final double currentMortgageRate = householdTerminal.getCurrentMortgageRate();
        final double futureAsset = xOpt.get(0);
        final double consumption = householdTerminal.getConsumption();
        final double wage = householdTerminal.getWage();
        final double newQuantity = 0;
        final double newMortgageRate = Double.NaN;
        final double newOriginationLTV = Double.NaN;
        final double status = householdTerminal.getStatus();
        final double averageBequest = (householdTerminal.getBequestUp() + householdTerminal.getBequestDown()) / 2;

        double[] output = {age, value, currentAsset, currentQuantitiy, currentMortgageRate, currentOriginationLTV, ficoScore, currentLTV,
                           wage, consumption, futureAsset, newQuantity, newMortgageRate, newOriginationLTV, status, averageBequest};
        return output;
    }

    double[] interim(int age, double currentAsset, double currentQuantity, double currentOriginationLTV,
                     double ficoScore, double currentLTV, double currentProductivity, Spline2DInterpolate creditSurface,
                     double[] productivityGrid, Grid gridFuture, DataVector alphaFuture){

        ObjectiveFunctionInterim householdInterim = new ObjectiveFunctionInterim();
        householdInterim.setAge(age);
        householdInterim.setCurrentAsset(currentAsset);
        householdInterim.setCurrentQuantity(currentQuantity);
        householdInterim.setCurrentOriginationLTV(currentOriginationLTV);
        householdInterim.setFicoScore(ficoScore);
        householdInterim.setCurrentLTV(currentLTV);
        householdInterim.setCurrentProductivity(currentProductivity);
        householdInterim.setCreditSurface(creditSurface);
        householdInterim.setProductivityGrid(productivityGrid);
        householdInterim.setGridFuture(gridFuture);
        householdInterim.setAlphaFuture(alphaFuture);

        // dimension of the domain
        final long dim = householdInterim.getNumberOfParameters();
        // B-spline degree
        final long p = 3;
        // max number of grid points
        final long N = 60;
        // adaptivity of grid generation
        final double gamma = 0.95;

        // First, we define a grid with modified B-spline basis functions and an iterative grid generator, which can generate the grid adaptively.
        Grid grid = Grid.createModBsplineGrid(dim, p);
        OptIterativeGridGeneratorRitterNovak gridGen = new OptIterativeGridGeneratorRitterNovak(householdInterim, grid, N, gamma);

        // With the iterative grid generator, we generate adaptively a sparse grid.
        gridGen.generate();

        // Then, we hierarchize the function values to get hierarchical B-spline coefficients of the B-spline sparse grid interpolant \hat{f} : [0, 1]^d -> R
        final DataVector functionValues = gridGen.getFunctionValues();
        DataVector coeffs = new DataVector(functionValues.getSize());
        HierarchisationSLE hierarchisationSLE = new HierarchisationSLE(grid);
        AutoSLESolver sleSolver = new AutoSLESolver();

        sleSolver.solve(hierarchisationSLE, gridGen.getFunctionValues(), coeffs);

        // Now define the interpolant $\hat{f}$ and its gradient, $\nabla \hat{f}$ for use with the gradient method (steepest descent).
        InterpolantScalarFunction householdInterimHat = new InterpolantScalarFunction(grid, coeffs);
        InterpolantScalarFunctionGradient householdInterimHatGradient = new InterpolantScalarFunctionGradient(grid, coeffs);
        OptGradientDescent gradientDescent = new OptGradientDescent(householdInterimHat, householdInterimHatGradient);
        DataVector x0 = new DataVector(dim);
        double fX0;

        GridStorage gridStorage = grid.getStorage();
        // index of grid point with minimal function value
        int x0Index = 0;
        fX0 = functionValues.get(0);
        for (int i = 1; i < functionValues.getSize(); i++) {
            if (functionValues.get(i) < fX0) {
                fX0 = functionValues.get(i);
                x0Index = i;
            }
        }
        x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));

        // We apply the gradient method and print the results.
        gradientDescent.setStartingPoint(x0);
        gradientDescent.optimize();
        DataVector xOpt = gradientDescent.getOptimalPoint();
        final double value = - householdInterim.eval(xOpt);
        final double futureAsset = householdInterim.getFutureAsset();
        final double consumption = householdInterim.getConsumption();
        final double currentMortgageRate = householdInterim.getCurrentMortgageRate();
        final double wage = householdInterim.getWage();
        final double newQuantity = householdInterim.getNewQuantity();
        final double newMortgageRate = householdInterim.getNewMortgageRate();
        final double newOriginationLTV = householdInterim.getNewOriginationLTV();
        final double status = householdInterim.getStatus();
        final double averageBequest = Double.NaN;

        double[] output = {age, value, currentAsset, currentQuantity, currentMortgageRate, currentOriginationLTV, ficoScore, currentLTV,
                           wage, consumption, futureAsset, newQuantity, newMortgageRate, newOriginationLTV, status, averageBequest};
        return output;
    }

    double[] newborn(int age, double currentAsset, double ficoScore, double currentProductivity, Spline2DInterpolate creditSurface,
                     double[] productivityGrid, Grid gridFuture, DataVector alphaFuture){

        ObjectiveFunctionNewborn householdNewborn = new ObjectiveFunctionNewborn();
        householdNewborn.setAge(age);
        householdNewborn.setCurrentAsset(currentAsset);
        householdNewborn.setFicoScore(ficoScore);
        householdNewborn.setCurrentProductivity(currentProductivity);
        householdNewborn.setCreditSurface(creditSurface);
        householdNewborn.setProductivityGrid(productivityGrid);
        householdNewborn.setGridFuture(gridFuture);
        householdNewborn.setAlphaFuture(alphaFuture);

        // dimension of the domain
        final long dim = householdNewborn.getNumberOfParameters();
        // B-spline degree
        final long p = 3;
        // max number of grid points
        final long N = 60;
        // adaptivity of grid generation
        final double gamma = 0.95;

        // First, define a grid with modified B-spline basis function and an iterative grid generator, which can generate the grid adaptively
        Grid grid = Grid.createModBsplineGrid(dim, p);
        OptIterativeGridGeneratorRitterNovak gridGen = new OptIterativeGridGeneratorRitterNovak(householdNewborn, grid, N, gamma);

        // With the iterative grid generator, we generate adaptively a sparse grid.
        gridGen.generate();

        // Now hierarchize to get coefficients of the B-spline grid
        final DataVector functionValues = gridGen.getFunctionValues();
        DataVector coeffs = new DataVector(functionValues.getSize());
        HierarchisationSLE hierarchisationSLE = new HierarchisationSLE(grid);
        AutoSLESolver sleSolver = new AutoSLESolver();

        sleSolver.solve(hierarchisationSLE, gridGen.getFunctionValues(), coeffs);
        // Now define the interpolant $\hat{f}$ and its gradient, $\nabla \hat{f}$ for use with the gradient method (steepest descent).
        InterpolantScalarFunction householdNewbornHat = new InterpolantScalarFunction(grid, coeffs);
        InterpolantScalarFunctionGradient householdNewbornHatGradient = new InterpolantScalarFunctionGradient(grid, coeffs);
        OptGradientDescent gradientDescent = new OptGradientDescent(householdNewbornHat, householdNewbornHatGradient);
        DataVector x0 = new DataVector(dim);
        double fX0;

        GridStorage gridStorage = grid.getStorage();
        // index of grid point with minimal function value
        int x0Index = 0;
        fX0 = functionValues.get(0);
        for (int i = 1; i < functionValues.getSize(); i++) {
            if (functionValues.get(i) < fX0) {
                fX0 = functionValues.get(i);
                x0Index = i;
            }
        }
        x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));

        // apply the gradient method and collect the result
        gradientDescent.setStartingPoint(x0);
        gradientDescent.optimize();
        DataVector xOpt = gradientDescent.getOptimalPoint();
        final double value = - householdNewborn.eval(xOpt);
        final double futureAsset = householdNewborn.getFutureAsset();
        final double consumption = householdNewborn.getConsumption();
        final double wage = householdNewborn.getWage();
        final double newQuantity = householdNewborn.getNewQuantity();
        final double newMortgageRate = householdNewborn.getNewMortgageRate();
        final double newOriginationLTV = householdNewborn.getNewOriginationLTV();
        final double status = householdNewborn.getStatus();

        // not applicable variables for newborns
        final double currentQuantity = 0;
        final double currentMortgageRate = Double.NaN;
        final double currentOriginationLTV = Double.NaN;
        final double currentLTV = Double.NaN;
        final double averageBequest = Double.NaN;

        double[] output = {age, value, currentAsset, currentQuantity, currentMortgageRate, currentOriginationLTV, ficoScore, currentLTV,
                wage, consumption, futureAsset, newQuantity, newMortgageRate, newOriginationLTV, status, averageBequest};
        return output;

    }

    public static void main(String[] args){

        sgpp.LoadJSGPPLib.loadJSGPPLib();
        sgpp.jsgpp.omp_set_num_threads(1);
        // verbosity of the output
        sgpp.Printer.getInstance().disableStatusPrinting();

        Control control = new Control();

        String creditSurfaceData = "/home/ahyan/Dropbox/Housing Market and Leverage Cycle/Code/src/Housing/Data/CreditSurface.csv";
        int numberOfFicos = 7, numberOfOltvs = 24;
        Spline2DInterpolate creditSurface = control.creditSurface(creditSurfaceData, numberOfFicos, numberOfOltvs);

        int ne = 9;
        double sigma_eps = 0.2058, lambda_eps = 0.99;
        double[] productivityGrid = new Tauchen().productivityGrid(ne, sigma_eps, lambda_eps, 3);

        //================TERMINAL: START================//

        int dimension = 6;
        Grid gridTerminal = Grid.createLinearGrid(dimension);
        GridStorage gridStorageTerminal = gridTerminal.getStorage();

        int level = 4;
        gridTerminal.getGenerator().regular(level);

        DataVector alphaTerminal = new DataVector(gridStorageTerminal.getSize());
        alphaTerminal.setAll(0.0);

        int Z = (int) gridStorageTerminal.getSize();
        System.out.println("Grid points: " + Z);

        int age = 30;
        System.out.println("age: " + age);

        ArrayList<double[]> outputMatrix = new ArrayList<>();

        for (int i = 0; i < gridStorageTerminal.getSize(); i++) {
            GridPoint gp = gridStorageTerminal.getPoint(i);
            double[] outputTerminal = control.terminal(age, gp.getStandardCoordinate(0), gp.getStandardCoordinate(1),
                    gp.getStandardCoordinate(2), gp.getStandardCoordinate(3),
                    gp.getStandardCoordinate(4), gp.getStandardCoordinate(5),
                    creditSurface, productivityGrid);
            outputMatrix.add(outputTerminal);
            alphaTerminal.set(i, outputTerminal[1]);
        }
        System.out.println(alphaTerminal.l2Norm());

        jsgpp.createOperationHierarchisation(gridTerminal).doHierarchisation(alphaTerminal);
        System.out.println(alphaTerminal.l2Norm());

        //================TERMINAL: END==================//

        //================INTERIM: START================//

        Grid gridFuture = gridTerminal;
        DataVector alphaFuture = alphaTerminal;

        Grid gridCurrent = Grid.createLinearGrid(dimension);
        GridStorage gridStorageCurrent = gridCurrent.getStorage();

        gridCurrent.getGenerator().regular(level);

        DataVector alphaCurrent = new DataVector(gridStorageCurrent.getSize());
        alphaCurrent.setAll(0.0);

        while (age > 2){
            age--;
            System.out.println("age: " + age);

            for (int i = 0; i < gridStorageCurrent.getSize(); i++){
                GridPoint gp = gridStorageCurrent.getPoint(i);
                double[] outputInterim = control.interim(age, gp.getStandardCoordinate(0), gp.getStandardCoordinate(1),
                        gp.getStandardCoordinate(2), gp.getStandardCoordinate(3),
                        gp.getStandardCoordinate(4), gp.getStandardCoordinate(5),
                        creditSurface, productivityGrid, gridFuture, alphaFuture);
                outputMatrix.add(outputInterim);
                alphaCurrent.set(i, outputInterim[1]);
            }

            jsgpp.createOperationHierarchisation(gridCurrent).doHierarchisation(alphaCurrent);

            alphaFuture = alphaCurrent;
            gridFuture = gridCurrent;

        }

        //================INTERIM: END==================//

        //================NEWBORN: START================//

        dimension = 3;
        Grid gridNewborn = Grid.createLinearGrid(dimension);
        GridStorage gridStorageNewborn = gridNewborn.getStorage();

        level = 6;
        gridNewborn.getGenerator().regular(level);

        DataVector alphaNewborn = new DataVector(gridStorageNewborn.getSize());
        alphaNewborn.setAll(0.0);

        Z = (int) gridStorageNewborn.getSize();
        System.out.println("Grid points: " + Z);

        age--;
        System.out.println("age: " + age);

        for (int i = 0; i < gridStorageNewborn.getSize(); i++){
            GridPoint gp = gridStorageNewborn.getPoint(i);
            double[] outputNewborn = control.newborn(age, gp.getStandardCoordinate(0), gp.getStandardCoordinate(1),
                    gp.getStandardCoordinate(2), creditSurface, productivityGrid, gridFuture, alphaFuture);
            outputMatrix.add(outputNewborn);
            alphaNewborn.set(i, outputNewborn[1]);
        }

        jsgpp.createOperationHierarchisation(gridNewborn).doHierarchisation(alphaNewborn);

        //================NEWBORN: END==================//

        new Linspace().dataWriter(outputMatrix, "src/Sparse/OutputData.csv");

        System.out.println("=================================");

        double demand = 0;
        double supply = 0;
        for (int i = 0; i < outputMatrix.size(); i++){
            supply = supply + outputMatrix.get(i)[3];
            demand = demand + outputMatrix.get(i)[11];
        }
        System.out.println("supply: " + supply);
        System.out.println("demand: " + demand);

    }
}
