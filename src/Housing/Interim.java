package Housing;

import com.imsl.math.Spline2DInterpolate;
import sgpp.*;

import static java.lang.Math.*;
import static java.lang.Math.max;

public class Interim {

    double[] householdInterim(int age, double currentAsset, double currentQuantity, double currentOriginationLTV,
                              double ficoScore, double currentLTV, double currentProductivity,
                              Spline2DInterpolate creditSurface, double[] assetGrid, double[] productivityGrid,
                              double[] quantityGrid, double[] originationLTVGrid, Grid gridFuture, DataVector alphaFuture){

        Grids grids = new Grids();

        double minimumConsumption = 1e-5, minimumQuantity = 1e-5, minimumAsset = 1e-5;
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
        double[] pi = {1d/2, 1d/2};

        // existing mortgage terms
        double currentMortgagePayment = currentMortgageRate * currentOriginationLTV * housePrice * currentQuantity / (1 - pow(1 + currentMortgageRate, -term));
        if (currentLTV == 0){currentMortgagePayment = 0;}
        double futureBalance = (1 + currentMortgageRate) * (1.2 * currentLTV) * housePrice * currentQuantity - currentMortgagePayment;
        double ltvUp = futureBalance / (housePriceUp * currentQuantity);
        double ltvDown = futureBalance / (housePriceDown * currentQuantity);

        // new mortgage terms
        double newMortgageRate, newMortgagePayment, newOriginationLTV, newQuantity;

        // Tauchen
        int ne = 9;
        double sigma_eps = 0.02058, lambda_eps = 0.99;
        double[] transitionVector = new Tauchen().transitionVector(currentProductivity, productivityGrid, ne, sigma_eps, lambda_eps);

        double currentProductivityTransformed = (productivityGrid[productivityGrid.length - 1] - productivityGrid[0]) * currentProductivity + productivityGrid[0];
        double wage = exp(currentProductivityTransformed);

        // value function terms
        double utility, futureAsset, consumption;
        double VV, CC = 0, AA = 0, OO = 0, QQ = 0, status = 1;
        double VVC = -1e5, CCC = 0, AAC = 0;
        double VVR = -1e5, CCR = 0, AAR = 0, OOR = 0;
        double VVS = -1e5, CCS = 0, AAS = 0, OOS = 0, QQS = 0;
        double VVD = -1e5, CCD = 0, AAD = 0, OOD = 0, QQD = 0;
        double expectedC, expectedR, expectedS, expectedD;
        double[] expectedCArray, expectedRArray, expectedSArray, expectedDArray;
        OperationEval operationEval = jsgpp.createOperationEval(gridFuture);

        // sundry cost terms
        double closingCostFactor = 0.05; double maintenanceCostFactor = 0.1; double movingCostFactor = 0; double defaultPenalty = 1;
        double closingCost; double maintenanceCost, movingCost;
        double windfall; double newDownPayment;

        for (double asset : assetGrid){

            futureAsset = asset;
            System.out.println("futureAsset: " + futureAsset);
            //#####################################################//
            //--------------STAY CURRENT: START--------------------//
            //#####################################################//
            expectedCArray = grids.expectations(futureAsset, currentQuantity, currentOriginationLTV, ficoScore, ltvUp, ltvDown,
                                                        productivityGrid, transitionVector, alphaFuture, operationEval);
            expectedC = pi[0] * expectedCArray[0] + pi[1] * expectedCArray[1];
            maintenanceCost = maintenanceCostFactor * currentMortgagePayment;
            consumption = (1 + interestRate) * currentAsset + wage - futureAsset - currentMortgagePayment - maintenanceCost;

            if (consumption <= 0){consumption = minimumConsumption;}
            utility = (1 - alpha) * log(consumption) + alpha * log(currentQuantity) + beta * expectedC;

            if (utility >= VVC){
                VVC = utility;
                CCC = consumption;
                AAC = futureAsset;
            }
            //#####################################################//
            //--------------STAY CURRENT: END----------------------//
            //#####################################################//

            for (double originationLTV : originationLTVGrid){
                //#################################################//
                //----------------REFINANCE: START-----------------//
                //#################################################//
                newOriginationLTV = originationLTV;
                newMortgageRate = creditSurface.value(ficoScore, newOriginationLTV);
                newMortgagePayment = newMortgageRate * newOriginationLTV * housePrice * currentQuantity / (1 - pow(1 + newMortgageRate, -term));
                windfall = (newOriginationLTV - currentLTV) * housePrice * currentQuantity;
                maintenanceCost = maintenanceCostFactor * newMortgagePayment;
                closingCost = closingCostFactor * housePrice * currentQuantity;
                futureBalance = (1 + newMortgageRate) * newOriginationLTV * currentQuantity * housePrice - newMortgagePayment;
                ltvUp = futureBalance / (housePriceUp * currentQuantity);
                ltvDown = futureBalance / (housePriceDown * currentQuantity);

                expectedRArray = grids.expectations(futureAsset, currentQuantity, newOriginationLTV, ficoScore, ltvUp, ltvDown,
                                                   productivityGrid, transitionVector, alphaFuture, operationEval);
                expectedR = pi[0] * expectedRArray[0] + pi[1] * expectedRArray[1];
                consumption = (1 + interestRate) * currentAsset + wage + windfall - futureAsset - newMortgagePayment - maintenanceCost - closingCost;

                if (consumption <= 0){consumption = minimumConsumption;}
                utility = (1 - alpha) * log(consumption) + alpha * log(currentQuantity) + beta * expectedR;

                if (utility >= VVR){
                    VVR = utility;
                    CCR = consumption;
                    AAR = futureAsset;
                    OOR = newOriginationLTV;
                }
                //##################################################//
                //--------------------REFINANCE: END----------------//
                //##################################################//

                for (double quantity : quantityGrid){
                    //##################################################//
                    //---------------BUY/SELL: START--------------------//
                    //##################################################//
                    newQuantity = quantity;
                    newMortgagePayment = newMortgageRate * newOriginationLTV * housePrice * newQuantity / (1 - pow(1 + newMortgageRate, -term));
                    newDownPayment = (1 - newOriginationLTV) * housePrice * newQuantity;
                    windfall = (1 - currentLTV) * housePrice * currentQuantity;
                    closingCost = closingCostFactor * housePrice * newQuantity;
                    maintenanceCost = maintenanceCostFactor * newMortgagePayment;
                    movingCost = movingCostFactor * newMortgagePayment;
                    futureBalance = (1 + newMortgageRate) * newOriginationLTV * newQuantity * housePrice - newMortgagePayment;
                    ltvUp = futureBalance / (housePriceUp * newQuantity);
                    ltvDown = futureBalance / (housePriceDown * newQuantity);

                    expectedSArray = grids.expectations(futureAsset, newQuantity, newOriginationLTV, ficoScore, ltvUp, ltvDown,
                                                        productivityGrid, transitionVector, alphaFuture, operationEval);
                    expectedS = pi[0] * expectedSArray[0] + pi[1] * expectedSArray[1];

                    consumption = (1 + interestRate) * currentAsset + wage + windfall
                                - futureAsset - newMortgagePayment - maintenanceCost - closingCost - movingCost - newDownPayment;

                    if (consumption <= 0){consumption = minimumConsumption;}
                    utility = (1 - alpha) * log(consumption) + alpha * log(newQuantity) + beta * expectedS;

                    if (utility >= VVS){
                        VVS = utility;
                        CCS = consumption;
                        AAS = futureAsset;
                        QQS = newQuantity;
                        OOS = newOriginationLTV;
                    }
                    //#####################################################//
                    //--------------------BUY/SELL: END--------------------//
                    //#####################################################//

                    //#####################################################//
                    //--------------------DEFAULT: START-------------------//
                    //#####################################################//
                    expectedDArray = grids.expectations(futureAsset, newQuantity, newOriginationLTV, 0, ltvUp, ltvDown,
                                                        productivityGrid, transitionVector, alphaFuture, operationEval);
                    expectedD = pi[0] * expectedDArray[0] + pi[1] * expectedDArray[1];
                    newDownPayment = (1 - newOriginationLTV) * housePrice * newQuantity;
                    consumption = (1 + interestRate) * currentAsset + wage - futureAsset - newMortgagePayment - newDownPayment - closingCost - movingCost;

                    if (consumption <= 0){consumption = minimumConsumption;}
                    utility = (1 - alpha) * log(consumption) + alpha * log(newQuantity) - defaultPenalty + beta * expectedD;

                    if (utility >= VVD){
                        VVD = utility;
                        CCD = consumption;
                        AAD = futureAsset;
                        QQD = newQuantity;
                        OOD = newOriginationLTV;
                    }
                    //########################################################//
                    //--------------------DEFAULT: END------------------------//
                    //########################################################//
                }
            }
        }

        System.out.println(VVC + " , " + VVR + " , " + VVS + " , " + VVD);
        VV = max(max(VVC, VVR), max(VVS, VVD));
        if (VV == VVC){
            CC = CCC;
            AA = AAC;
            QQ = currentQuantity;
            OO = currentOriginationLTV;
            status = 1;
        }else if (VV == VVR){
            CC = CCR;
            AA = AAR;
            QQ = currentQuantity;
            OO = OOR;
            status = 2;
        }else if (VV == VVS){
            CC = CCS;
            AA = AAS;
            QQ = QQS;
            OO = OOS;
            status = 3;
        }else if (VV == VVD){
            CC = CCD;
            AA = AAD;
            QQ = QQD;
            OO = OOD;
            status = 4;
        }

        /*
         * header = { "age", "value", "current_asset", "current_quantity",
         *            "current_origination_ltv", "fico_score", "current_ltv",
         *            "wage", "consumption", "future_asset", "future_quantity",
         *            "future_oltv", "status" };
         */
        double[] output = {age, VV, currentAsset, currentQuantity, currentOriginationLTV,
                ficoScore, currentLTV, wage, CC, AA, QQ, OO, status};
        return output;

    }
}


