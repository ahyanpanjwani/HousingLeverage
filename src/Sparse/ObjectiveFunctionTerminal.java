package Sparse;

import com.imsl.math.Spline2DInterpolate;

import static java.lang.Math.*;

public class ObjectiveFunctionTerminal extends sgpp.ScalarFunction {

    public ObjectiveFunctionTerminal(){super(1);}

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

    double housePrice;

    // outputs
    double currentMortgageRate;
    double consumption;
    double wage;
    double status = 5.0;
    double bequestUp;
    double bequestDown;


    public void setAge(int age) { this.age = age; }

    public void setCurrentAsset(double currentAsset){
        this.currentAsset = currentAsset;
    }

    public void setCurrentQuantity(double currentQuantity){ this.currentQuantity = currentQuantity; }

    public void setCurrentOriginationLTV(double currentOriginationLTV){ this.currentOriginationLTV = currentOriginationLTV;}

    public void setFicoScore(double ficoScore) { this.ficoScore = ficoScore; }

    public void setCurrentLTV(double currentLTV) { this.currentLTV = currentLTV; }

    public void setCurrentProductivity(double currentProductivity) { this.currentProductivity = currentProductivity; }

    public void setCreditSurface(Spline2DInterpolate creditSurface) { this.creditSurface = creditSurface; }

    public void setProductivityGrid(double[] productivityGrid) { this.productivityGrid = productivityGrid; }

    public double eval(sgpp.DataVector x){

        double futureAsset = 1 * x.get(0);

        double minimumConsumption = 0.0001;
        double minimumQuantity = 0.0001;
        double minimumBequest = 0.0001;
        if (currentQuantity == 0){currentQuantity = minimumQuantity;}

        // set params
        double alpha = 0.33, beta = 0.97;
        housePrice = 40;
        double houseDrift = 0.0, houseVol = 0.06;
        double housePriceUp = housePrice * exp(houseDrift + houseVol);
        double housePriceDown = housePrice * exp(houseDrift - houseVol);
        double interestRate = 0.03;
        double term = 30;
        currentMortgageRate = creditSurface.value(ficoScore, currentOriginationLTV);

        // probabilities of the two states, 1/2
        double piUp = 0.5, piDown = 1 - piUp;

        // mortgage terms
        double currentMortgagePayment = currentMortgageRate * currentOriginationLTV * housePrice * currentQuantity / (1 - pow(1 + currentMortgageRate, -term));
        if (currentLTV == 0){currentMortgagePayment = 0;}
        double futureBalance = (1 + currentMortgageRate) * (currentLTV) * housePrice * currentQuantity - currentMortgagePayment;
        double ltvUp = futureBalance / (housePriceUp * currentQuantity);
        double ltvDown = futureBalance / (housePriceDown * currentQuantity);

        // value function terms
        double utility, expected;

        // productivity
        double[] productivityGrid = {-4.376632559721463, -3.282474419791097, -2.1883162798607314, -1.0941581399303657, 0.0, 1.0941581399303657, 2.1883162798607314, 3.282474419791097, 4.376632559721463};
        double currentProductivityTransformed = (productivityGrid[productivityGrid.length - 1] - productivityGrid[0]) * currentProductivity + productivityGrid[0];
        wage = exp(currentProductivityTransformed);

        bequestUp = (1 + interestRate) * futureAsset + (1 - ltvUp) * housePriceUp * currentQuantity;
        bequestDown = (1 + interestRate) * futureAsset + (1 - ltvDown) * housePriceDown * currentQuantity;
        if (bequestUp <= 0){bequestUp = minimumBequest;}
        if (bequestDown <= 0){bequestDown = minimumBequest;}

        expected = piUp * log(bequestUp) + piDown * log(bequestDown);
        consumption = (1 + interestRate) * currentAsset + wage - currentMortgagePayment - futureAsset;
        if (consumption <= 0){consumption = minimumConsumption;}

        utility = (1 - alpha) * log(consumption) + alpha * log(currentQuantity) + beta * expected;
        //System.out.println("u:" + utility);
        return -utility;
    }

    public double getConsumption() {
        return consumption;
    }

    public double getCurrentMortgageRate() { return currentMortgageRate; }

    public double getWage() { return wage; }

    public double getStatus() { return  status; }

    public double getBequestUp() { return bequestUp; }

    public double getBequestDown() {return  bequestDown; }

    //public void clone(sgpp.SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunction_t clone) {}

    public void clone(sgpp.SWIGTYPE_p_std__unique_ptrT_sgpp__optimization__optimizer__UnconstrainedOptimizer_t clone) {}

}
