package Housing;

import com.imsl.math.Spline2DInterpolate;

import static java.lang.Math.*;

public class Terminal {

    double[] householdTerminal(int age, double currentAsset, double currentQuantity, double currentOriginationLTV,
                               double ficoScore, double currentLTV, double currentProductivity,
                               Spline2DInterpolate creditSurface, double[] assetGrid, double[] productivityGrid){

        double minimumConsumption = 0.001, minimumQuantity = 0.001, minimumAsset = 0.001;
        if (currentAsset == 0){currentAsset = minimumAsset;}
        if (currentQuantity == 0){currentQuantity = minimumQuantity;}

        // set params
        double alpha = 0.33, beta = 0.97;
        double housePrice = 50, houseDrift = 0.0, houseVol = 0.06;
        double housePriceUp = housePrice * exp(houseDrift + houseVol);
        double housePriceDown = housePrice * exp(houseDrift - houseVol);
        double interestRate = 0.03;
        double term = 30;
        double currentMortgageRate = creditSurface.value(ficoScore, currentOriginationLTV);

        // probabilities of the two states, 1/2
        double piUp = 0.5, piDown = 1 - piUp;

        // mortgage terms
        double currentMortgagePayment = currentMortgageRate * currentOriginationLTV * housePrice * currentQuantity / (1 - pow(1 + currentMortgageRate, -term));
        if (currentLTV == 0){currentMortgagePayment = 0;}
        double futureBalance = (1 + currentMortgageRate) * (1.2 * currentLTV) * housePrice * currentQuantity - currentMortgagePayment;
        double ltvUp = futureBalance / (housePriceUp * currentQuantity);
        double ltvDown = futureBalance / (housePriceDown * currentQuantity);

        // value function terms
        double utility, futureAsset, consumption, expected;
        double bequestUp, bequestDown;
        double VV = -10;
        double CC = 0, QQ = 0, AA = 0, status = 5;

        // productivity
        double currentProductivityTransformed = (productivityGrid[productivityGrid.length - 1] - productivityGrid[0]) * currentProductivity + productivityGrid[0];
        double wage = exp(currentProductivityTransformed);

        for (double asset : assetGrid) {

            futureAsset = 4 * asset;
            bequestUp = (1 + interestRate) * futureAsset + (1 - ltvUp) * housePriceUp * currentQuantity;
            bequestDown = (1 + interestRate) * futureAsset + (1 - ltvDown) * housePriceDown * currentQuantity;

            expected = piUp * log(bequestUp) + piDown * log(bequestDown);
            consumption = (1 + interestRate) * currentAsset + wage - currentMortgagePayment - futureAsset;
            if (consumption <= 0){consumption = minimumConsumption;}

            // calculate the utility: flow utility (log/Cobb-Douglas) plus the expected log of bequest
            utility = (1 - alpha) * log(consumption) + alpha * log(currentQuantity) + beta * expected;

            if (utility >= VV) {
                VV = utility;
                CC = consumption;
                QQ = currentQuantity;
                AA = futureAsset / 4;
            }

        }

        /**
         * header = { "age", "value", "current_asset", "current_quantity",
         *            "current_origination_ltv", "fico_score", "current_ltv",
         *            "wage", "consumption", "future_asset", "future_quantity",
         *            "future_oltv", "status" };
         */
        double[] output = {age, VV, currentAsset, currentQuantity, currentOriginationLTV,
                           ficoScore, currentLTV, wage, CC, AA, QQ, currentOriginationLTV, status};
        return output;
    }

}
