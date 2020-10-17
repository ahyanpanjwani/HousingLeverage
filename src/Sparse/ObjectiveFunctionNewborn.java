package Sparse;

import com.imsl.math.Spline2DInterpolate;
import sgpp.*;

import static java.lang.Math.*;

public class ObjectiveFunctionNewborn extends ScalarFunction {

    // inputs
    int age;
    double currentAsset;
    double ficoScore;
    double currentProductivity;
    Spline2DInterpolate creditSurface;
    double[] productivityGrid;
    Grid gridFuture;
    DataVector alphaFuture;

    // outputs
    double utility;
    double futureAsset;
    double newMortgageRate;
    double newOriginationLTV;
    double newQuantity;
    double wage;

    double consumption;
    double status;

    Tauchen tauchen = new Tauchen();
    OperationEval operationEval;

    public void setAge(int age) {
        this.age = age;
    }

    public void setCurrentAsset(double currentAsset) {
        this.currentAsset = currentAsset;
    }

    public void setFicoScore(double ficoScore) {
        this.ficoScore = ficoScore;
    }

    public void setCurrentProductivity(double currentProductivity) {
        this.currentProductivity = currentProductivity;
    }

    public void setCreditSurface(Spline2DInterpolate creditSurface) {
        this.creditSurface = creditSurface;
    }

    public void setProductivityGrid(double[] productivityGrid) {
        this.productivityGrid = productivityGrid;
    }

    public void setGridFuture(Grid gridFuture) {
        this.gridFuture = gridFuture;
    }

    public void setAlphaFuture(DataVector alphaFuture) {
        this.alphaFuture = alphaFuture;
    }

    public ObjectiveFunctionNewborn(){super(3);}

    public double eval(DataVector x){

        // the choice variables
        futureAsset = x.get(0);
        newOriginationLTV = x.get(1);
        newQuantity = x.get(2);

        operationEval = jsgpp.createOperationEval(gridFuture);

        // enforce bounds to avoid nans
        double minimumConsumption = 0.0001;
        double minimumQuantity = 0.0001;
        if (newQuantity == 0){newQuantity = minimumQuantity;}

        // set params
        double alpha = 0.33, beta = 0.97;
        double housePrice = 40, houseDrift = 0.0, houseVol = 0.06;
        double housePriceUp = housePrice * exp(houseDrift + houseVol);
        double housePriceDown = housePrice * exp(houseDrift - houseVol);
        double interestRate = 0.03;
        double term = 30;
        newMortgageRate = creditSurface.value(ficoScore, newOriginationLTV);

        // Tauchen
        int ne = 9;
        double sigma_eps = 0.2058, lambda_eps = 0.99;
        double[] transitionVector = tauchen.transitionVector(currentProductivity, productivityGrid, ne, sigma_eps, lambda_eps);

        double currentProductivityTransformed = (productivityGrid[productivityGrid.length - 1] - productivityGrid[0]) * currentProductivity + productivityGrid[0];
        wage = exp(currentProductivityTransformed);

        // sundry cost terms
        double closingCostFactor = 0.05; double maintenanceCostFactor = 0.1; double movingCostFactor = 0.25;
        double closingCost; double maintenanceCost, movingCost;
        double windfall; double newDownPayment; double newMortgagePayment;

        // probabilities of the two states, 1/2
        double[] pi = {0.5, 0.5};

        newMortgagePayment = newMortgageRate * newOriginationLTV * housePrice * newQuantity / (1 - pow(1 + newMortgageRate, -term));
        newDownPayment = (1 - newOriginationLTV) * housePrice * newQuantity;
        closingCost = closingCostFactor * housePrice * newQuantity;
        maintenanceCost = maintenanceCostFactor * newMortgagePayment;
        movingCost = movingCostFactor * newMortgagePayment;

        double futureBalance = (1 + newMortgageRate) * newOriginationLTV * newQuantity * housePrice - newMortgagePayment;
        double ltvUp = futureBalance / (housePriceUp * newQuantity);
        double ltvDown = futureBalance / (housePriceDown * newQuantity);

        double[] expectedArray = {0, 0}; double expected;
        DataVector pUp = new DataVector(6);
        DataVector pDown = new DataVector(6);


        pUp.set(0, futureAsset); pUp.set(1, newQuantity); pUp.set(2, newOriginationLTV);
        pUp.set(3, ficoScore); pUp.set(4, ltvUp);
        pDown.set(0, futureAsset); pDown.set(1, newQuantity); pDown.set(2, newOriginationLTV);
        pDown.set(3, ficoScore); pDown.set(4, ltvDown);

        for (int iep = 0; iep < productivityGrid.length; iep++){
            double productivity = iep / (double) (productivityGrid.length - 1);
            pUp.set(5, productivity);
            pDown.set(5, productivity);
            expectedArray[0] = expectedArray[0] + transitionVector[iep] * operationEval.eval(alphaFuture, pUp);
            expectedArray[1] = expectedArray[1] + transitionVector[iep] * operationEval.eval(alphaFuture, pDown);
        }

        expected = pi[0] * expectedArray[0] + pi[1] * expectedArray[1];
        consumption = currentAsset + wage - newMortgagePayment - maintenanceCost - closingCost
                    - movingCost - newDownPayment - futureAsset;
        if (consumption <= 0){consumption = minimumConsumption;}

        utility = (1 - alpha) * log(consumption) + alpha * log(newQuantity) + beta * expected;

        return  - utility;
    }

    public void clone(sgpp.SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunction_t clone) {
    }

    public double getFutureAsset() {
        return futureAsset;
    }

    public double getNewMortgageRate() {
        return newMortgageRate;
    }

    public double getConsumption() {
        return consumption;
    }

    public double getWage() {
        return wage;
    }

    public double getNewOriginationLTV() {
        return newOriginationLTV;
    }

    public double getNewQuantity() {
        return newQuantity;
    }

    public double getStatus() {
        return status;
    }

}
