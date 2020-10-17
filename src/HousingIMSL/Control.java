package HousingIMSL;

import com.imsl.math.NelderMead;
import com.imsl.math.Spline2DInterpolate;
import sgpp.*;

public class Control {

    Spline2DInterpolate creditSurface(String filePathForCS, int nFICO, int nOLTV){
        CreditSurfaceInterpolation creditSurfaceInterpolation = new CreditSurfaceInterpolation();
        double[] ficos = creditSurfaceInterpolation.ficos();
        double[] originationLTVs = creditSurfaceInterpolation.originationLTVs(filePathForCS, nOLTV);
        double[][] surface = creditSurfaceInterpolation.creditSurface(filePathForCS, nFICO, nOLTV);
        Spline2DInterpolate interpolate = new Spline2DInterpolate(ficos, originationLTVs, surface, 2, 2);
        return interpolate;
    }

    double householdTerminal(double currentAsset, double quantity, double originationLTV, double fico, double currentLTV,
                             double productivity, Spline2DInterpolate creditSurface){
        Terminal terminal = new Terminal();
        terminal.setCurrentAsset(currentAsset);
        terminal.setQuantity(quantity);
        terminal.setOriginationLTV(originationLTV);
        terminal.setFico(fico);
        terminal.setCurrentLTV(currentLTV);
        terminal.setProductivity(productivity);
        terminal.setCreditSurface(creditSurface);

        NelderMead.Function objective = terminal.objectiveTerminal;
        NelderMead nelderMead = new NelderMead(objective, 1);
        nelderMead.setTolerance(1e-3);
        //double[] guess = {0.5};
        //nelderMead.setGuess(guess);
        //Random seed = new Random(0);
        //nelderMead.setRandomObject(seed);
        nelderMead.solve();

        return - nelderMead.getObjectiveValue();
    }

    public static void main(String[] args){
        Control control = new Control();
        Spline2DInterpolate creditSurface = control.creditSurface("src/HousingIMSL/CreditSurface.csv", 7, 24);

        Terminal terminal = new Terminal();
        double value = control.householdTerminal(0.5, 0.9, 0.97, 0, 0.9, 0.5,
                creditSurface);
        System.out.println("Direct: " + value);

        sgpp.LoadJSGPPLib.loadJSGPPLib();

        int dimension = 6;
        Grid grid = Grid.createLinearBoundaryGrid(dimension);
        GridStorage gridStorage = grid.getStorage();

        int level = 2;
        grid.getGenerator().regular(level);
        System.out.println("Number of points (CSG): " + gridStorage.getSize());

        DataVector alpha = new DataVector(gridStorage.getSize());
        alpha.setAll(0.0);

        for (int step = 0; step < 1; step++){
            for (int i = 0; i < gridStorage.getSize(); i++) {
                GridPoint gp = gridStorage.getPoint(i);
                alpha.set(i, control.householdTerminal(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1),
                        gp.getStandardCoordinate(2), gp.getStandardCoordinate(3), gp.getStandardCoordinate(4),
                        gp.getStandardCoordinate(5), creditSurface));
            }
            jsgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);
            SurplusRefinementFunctor functor = new SurplusRefinementFunctor(alpha, 1);
            grid.getGenerator().refine(functor);
            System.out.println("refinement step " + (step + 1) + ", new grid size: " + alpha.getSize());
            alpha.resizeZero(gridStorage.getSize());

        }

        OperationEval operationEval = jsgpp.createOperationEval(grid);

        DataVector p = new DataVector(6);
        p.set(0, 0.5);
        p.set(1, 0.9);
        p.set(2, 0.97);
        p.set(3, 0.0);
        p.set(4, 0.9);
        p.set(5, 0.5);
        System.out.println(operationEval.eval(alpha, p));

    }


}
