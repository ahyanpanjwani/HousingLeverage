package Sparse;

import com.imsl.math.Spline2DInterpolate;
import sgpp.*;

import static java.lang.Math.*;

public class ObjectiveFunctionInterim extends ScalarFunction {

    // inputs
    int age;
    double currentAsset;
    double currentQuantity;
    double currentOriginationLTV;
    double ficoScore;
    double currentLTV;
    double currentProductivity;
    Spline2DInterpolate creditSurface;
    double[] productivityGrid;
    Grid gridFuture;
    DataVector alphaFuture;

    // outputs
    double utility;
    double futureAsset;
    double currentMortgageRate;
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

    public void setCurrentQuantity(double currentQuantity) {
        this.currentQuantity = currentQuantity;
    }

    public void setCurrentOriginationLTV(double currentOriginationLTV) {
        this.currentOriginationLTV = currentOriginationLTV;
    }

    public void setFicoScore(double ficoScore) {
        this.ficoScore = ficoScore;
    }

    public void setCurrentLTV(double currentLTV) {
        this.currentLTV = currentLTV;
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

    public ObjectiveFunctionInterim(){super(3);}

    public double eval(sgpp.DataVector x){

        // the choice variables
        futureAsset = x.get(0);
        newOriginationLTV = x.get(1);
        newQuantity = x.get(2);

        operationEval = jsgpp.createOperationEval(gridFuture);

        // enforce bounds to avoid nans
        double minimumConsumption = 0.0001;
        double minimumQuantity = 0.0001;
        if (currentQuantity == 0){currentQuantity = minimumQuantity;}
        if (newQuantity == 0){newQuantity = minimumQuantity;}

        // set params
        double alpha = 0.33, beta = 0.97;
        double housePrice = 40, houseDrift = 0.0, houseVol = 0.06;
        double housePriceUp = housePrice * exp(houseDrift + houseVol);
        double housePriceDown = housePrice * exp(houseDrift - houseVol);
        double interestRate = 0.03;
        double term = 30;
        currentMortgageRate = creditSurface.value(ficoScore, currentOriginationLTV);

        // probabilities of the two states, 1/2
        double[] pi = {0.5, 0.5};

        // mortgage terms
        double currentMortgagePayment = currentMortgageRate * currentOriginationLTV * housePrice * currentQuantity / (1 - pow(1 + currentMortgageRate, -term));
        if (currentLTV == 0){currentMortgagePayment = 0;}
        double futureBalance = (1 + currentMortgageRate) * (currentLTV) * housePrice * currentQuantity - currentMortgagePayment;
        double ltvUp = futureBalance / (housePriceUp * currentQuantity);
        double ltvDown = futureBalance / (housePriceDown * currentQuantity);

        // Tauchen
        int ne = 9;
        double sigma_eps = 0.2058, lambda_eps = 0.99;
        double[] transitionVector = tauchen.transitionVector(currentProductivity, productivityGrid, ne, sigma_eps, lambda_eps);

        double currentProductivityTransformed = (productivityGrid[productivityGrid.length - 1] - productivityGrid[0]) * currentProductivity + productivityGrid[0];
        wage = exp(currentProductivityTransformed);

        // sundry cost terms
        double closingCostFactor = 0.05; double maintenanceCostFactor = 0.1; double movingCostFactor = 0.25;
        double closingCost; double maintenanceCost, movingCost;
        double windfall; double newDownPayment;
        double defaultPenalty = 1; double ficoScoreDefault = 0;

        // utility and consumption terms
        double expectedCurrent, expectedRefinance, expectedSell, expectedDefault;
        double utilityCurrent, utilityRefinance, utilitySell, utilityDefault;
        double consumptionCurrent, consumptionRefinance, consumptionSell, consumptionDefault;
        double[] expectedCurrentArray = {0, 0};
        double[] expectedRefinanceArray = {0, 0};
        double[] expectedSellArray = {0, 0};
        double[] expectedDefaultArray = {0, 0};
        DataVector pUp = new DataVector(6);
        DataVector pDown = new DataVector(6);

        //#####################################################//
        //--------------STAY CURRENT: START--------------------//
        //#####################################################//
        pUp.set(0, futureAsset); pUp.set(1, currentQuantity); pUp.set(2, currentOriginationLTV);
        pUp.set(3, ficoScore); pUp.set(4, ltvUp);
        pDown.set(0, futureAsset); pDown.set(1, currentQuantity); pDown.set(2, currentOriginationLTV);
        pDown.set(3, ficoScore); pDown.set(4, ltvDown);

        for (int iep = 0; iep < productivityGrid.length; iep++){
            double productivity = iep / (double) (productivityGrid.length - 1);
            pUp.set(5, productivity);
            pDown.set(5, productivity);
            expectedCurrentArray[0] = expectedCurrentArray[0] + transitionVector[iep] * operationEval.eval(alphaFuture, pUp);
            expectedCurrentArray[1] = expectedCurrentArray[1] + transitionVector[iep] * operationEval.eval(alphaFuture, pDown);
        }

        expectedCurrent = pi[0] * expectedCurrentArray[0] + pi[1] * expectedCurrentArray[1];
        maintenanceCost = maintenanceCostFactor * currentMortgagePayment;
        consumptionCurrent = (1 + interestRate) * currentAsset + wage - futureAsset - currentMortgagePayment - maintenanceCost;
        if (consumptionCurrent <= 0){consumptionCurrent = minimumConsumption;}
        utilityCurrent = (1 - alpha) * log(consumptionCurrent) + alpha * log(currentQuantity) + beta * expectedCurrent;
        //#####################################################//
        //--------------STAY CURRENT: END----------------------//
        //#####################################################//


        //#####################################################//
        //----------------REFINANCE: START---------------------//
        //#####################################################//
        newMortgageRate = creditSurface.value(ficoScore, newOriginationLTV);
        double newMortgagePayment = newMortgageRate * newOriginationLTV * housePrice * currentQuantity / (1 - pow(1 + newMortgageRate, -term));
        windfall = (newOriginationLTV - currentLTV) * housePrice * currentQuantity;
        maintenanceCost = maintenanceCostFactor * newMortgagePayment;
        closingCost = closingCostFactor * housePrice * currentQuantity;
        futureBalance = (1 + newMortgageRate) * newOriginationLTV * currentQuantity * housePrice - newMortgagePayment;
        ltvUp = futureBalance / (housePriceUp * currentQuantity);
        ltvDown = futureBalance / (housePriceDown * currentQuantity);

        pUp.set(0, futureAsset); pUp.set(1, currentQuantity); pUp.set(2, newOriginationLTV);
        pUp.set(3, ficoScore); pUp.set(4, ltvUp);
        pDown.set(0, futureAsset); pDown.set(1, currentQuantity); pDown.set(2, newOriginationLTV);
        pDown.set(3, ficoScore); pDown.set(4, ltvDown);

        for (int iep = 0; iep < productivityGrid.length; iep++){
            double productivity = iep / (double) (productivityGrid.length - 1);
            pUp.set(5, productivity);
            pDown.set(5, productivity);
            expectedRefinanceArray[0] = expectedRefinanceArray[0] + transitionVector[iep] * operationEval.eval(alphaFuture, pUp);
            expectedRefinanceArray[1] = expectedRefinanceArray[1] + transitionVector[iep] * operationEval.eval(alphaFuture, pDown);
        }

        expectedRefinance = pi[0] * expectedRefinanceArray[0] + pi[1] * expectedRefinanceArray[1];
        consumptionRefinance = (1 + interestRate) * currentAsset + wage + windfall - newMortgagePayment - maintenanceCost - closingCost;
        if (consumptionRefinance <= 0){consumptionRefinance = minimumConsumption;}
        utilityRefinance = (1 - alpha) * log(consumptionRefinance) + alpha * log(currentQuantity) + beta * expectedRefinance;
        //#####################################################//
        //----------------REFINANCE: END-----------------------//
        //#####################################################//


        //#####################################################//
        //----------------BUY/SELL: START----------------------//
        //#####################################################//
        newMortgagePayment = newMortgageRate * newOriginationLTV * housePrice * newQuantity / (1 - pow(1 + newMortgageRate, -term));
        newDownPayment = (1 - newOriginationLTV) * housePrice * newQuantity;
        windfall = (1 - currentLTV) * housePrice * currentQuantity;
        closingCost = closingCostFactor * housePrice * newQuantity;
        maintenanceCost = maintenanceCostFactor * newMortgagePayment;
        movingCost = movingCostFactor * newMortgagePayment;
        futureBalance = (1 + newMortgageRate) * newOriginationLTV * newQuantity * housePrice - newMortgagePayment;
        ltvUp = futureBalance / (housePriceUp * newQuantity);
        ltvDown = futureBalance / (housePriceDown * newQuantity);

        pUp.set(0, futureAsset); pUp.set(1, newQuantity); pUp.set(2, newOriginationLTV);
        pUp.set(3, ficoScore); pUp.set(4, ltvUp);
        pDown.set(0, futureAsset); pDown.set(1, newQuantity); pDown.set(2, newOriginationLTV);
        pDown.set(3, ficoScore); pDown.set(4, ltvDown);

        for (int iep = 0; iep < productivityGrid.length; iep++){
            double productivity = iep / (double) (productivityGrid.length - 1);
            pUp.set(5, productivity);
            pDown.set(5, productivity);
            expectedSellArray[0] = expectedSellArray[0] + transitionVector[iep] * operationEval.eval(alphaFuture, pUp);
            expectedSellArray[1] = expectedSellArray[1] + transitionVector[iep] * operationEval.eval(alphaFuture, pDown);
        }

        expectedSell = pi[0] * expectedSellArray[0] + pi[1] * expectedSellArray[1];
        consumptionSell = (1 + interestRate) * currentAsset + wage + windfall
                        - futureAsset - newMortgagePayment - maintenanceCost - closingCost - movingCost - newDownPayment;
        if (consumptionSell <= 0){consumptionSell = minimumConsumption;}
        utilitySell = (1 - alpha) * log(consumptionSell) + alpha * log(newQuantity) + beta * expectedSell;

        //#####################################################//
        //------------------BUY/SELL: END----------------------//
        //#####################################################//


        //#####################################################//
        //------------------DEFAULT: START---------------------//
        //#####################################################//
        newMortgageRate = creditSurface.value(ficoScoreDefault, newOriginationLTV);
        newMortgagePayment = newMortgageRate * newOriginationLTV * housePrice * newQuantity / (1 - pow(1 + newMortgageRate, -term));
        newDownPayment = (1 - newOriginationLTV) * housePrice * newQuantity;
        closingCost = closingCostFactor * housePrice * newQuantity;
        movingCost = movingCostFactor * newMortgagePayment;
        futureBalance = (1 + newMortgageRate) * newOriginationLTV * newQuantity * housePrice;
        ltvUp = futureBalance / (housePriceUp * newQuantity);
        ltvDown = futureBalance / (housePriceDown * newQuantity);

        pUp.set(0, futureAsset); pUp.set(1, newQuantity); pUp.set(2, newOriginationLTV);
        pUp.set(3, ficoScoreDefault); pUp.set(4, ltvUp);
        pDown.set(0, futureAsset); pDown.set(1, newQuantity); pDown.set(2, newOriginationLTV);
        pDown.set(3, ficoScoreDefault); pDown.set(4, ltvDown);

        for (int iep = 0; iep < productivityGrid.length; iep++){
            double productivity = iep / (double) (productivityGrid.length - 1);
            pUp.set(5, productivity);
            pDown.set(5, productivity);
            expectedDefaultArray[0] = expectedDefaultArray[0] + transitionVector[iep] * operationEval.eval(alphaFuture, pUp);
            expectedDefaultArray[1] = expectedDefaultArray[1] + transitionVector[iep] * operationEval.eval(alphaFuture, pDown);
        }

        expectedDefault = pi[0] * expectedDefaultArray[0] + pi[1] * expectedDefaultArray[1];

        consumptionDefault = (1 + interestRate) * currentAsset + wage - futureAsset - newDownPayment - closingCost - movingCost;
        if (consumptionDefault <= 0){consumptionDefault = minimumConsumption;}

        utilityDefault = (1 - alpha) * log(consumptionDefault) + alpha * log(newQuantity) - defaultPenalty + beta * expectedDefault;

        //#####################################################//
        //------------------DEFAULT: END-----------------------//
        //#####################################################//

        utility = max(max(utilityCurrent, utilityRefinance), max(utilitySell, utilityDefault));

        if (utility == utilityCurrent){
            consumption = consumptionCurrent;
            newMortgageRate = Double.NaN;
            newOriginationLTV = Double.NaN;
            newQuantity = currentQuantity;
            status = 1;
        }else if (utility == utilityRefinance){
            consumption = consumptionRefinance;
            newQuantity = currentQuantity;
            status = 2;
        }else if (utility == utilitySell){
            consumption = consumptionSell;
            status = 3;
        }else if (utility == utilityDefault){
            consumption = consumptionDefault;
            status = 4;
        }

       return - utility;
    }

    public void clone(sgpp.SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunction_t clone) {
    }

    public double getFutureAsset() { return futureAsset; }

    public double getCurrentMortgageRate() { return currentMortgageRate; }

    public double getConsumption() {
        return consumption;
    }

    public double getWage() { return wage; }

    public double getNewMortgageRate(){ return  newMortgageRate; }

    public double getNewOriginationLTV() { return newOriginationLTV; }

    public double getNewQuantity() { return newQuantity; }

    public double getStatus() {
        return status;
    }
}
