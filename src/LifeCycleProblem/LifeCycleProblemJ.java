package LifeCycleProblem;

import static java.lang.Math.*;

/**
 * This version of the model (v.J) will include FICO scores in the state space (so far we have only looked at FICO = 500)
 * in addition to having price uncertainty (developed in v.H.) and productivity uncertainty.
 *
 * Please note there is no v.I. in order to avoid typographic misunderstandings.
 *
 * This is the first iteration after improving the credit surface (Pegasus for root finding) and significant clean-up
 * using the Back-up package. Also it is the first update since we got cluster access.
 *
 * July 27, 2020
 */
public class LifeCycleProblemJ {

    double[] value_T(int age, int ia, int iq, int iOLTV, int iLTV, int iFICO, int ie, int T, int na,
                     double[] agrid, double[] qgrid, double[] oltvGrid, double[][] creditSurface,
                     double[][] ltvSchedule, double[] egrid, double alpha, double beta, double w,
                     double r, double pH, double muH, double sigmaH, int N){

        // mortgage terms
        double mortgageRate;
        double originationLTV;
        double quantity = 0;
        double mortgagePayment;
        double currentLTV;
        double balancePrime;       // balance in the next period after making this period's payment
        double ltvU, ltvD;

        // price expectations terms
        double delta_x = sqrt(pow(sigmaH, 2) + pow(muH, 2));      // step for house price
        double piU = 0.5 + (0.5 * muH / delta_x);                 // prob. of house price going up
        double piD = 0.5 - (0.5 * muH / delta_x);                 // prob. of house price going down
        double pHU = pH * exp(delta_x);                          // house price at the up node
        double pHD = pH * exp(-delta_x);                         // house price at the down node

        // value function terms
        double utility;
        double consumption;
        double bequestU;
        double bequestD;
        double expected;
        double VV = pow(-10, 5);
        double CC = 0;
        double QQ = 0;
        double AA = 0;

        if (age == T - 1){

            // set mortgage terms
            mortgageRate = creditSurface[iFICO][iOLTV];
            originationLTV = oltvGrid[iOLTV];
            quantity = qgrid[iq];
            mortgagePayment = mortgageRate * originationLTV * pH * quantity / (1 - pow(1 + mortgageRate, - N));
            currentLTV = ltvSchedule[iOLTV][iLTV];

            // if the agent own the house outright, then there is no mortgage payment to make
            if (currentLTV == 0){mortgagePayment = 0;}

            /* balance and ltv in the next period (and states) after making a periodic payment in the current period
             * (max ensures this is not negative; needed for guys who will pay off mortgage in current period)
             */
            balancePrime = (1 + mortgageRate) * currentLTV * pH * quantity - mortgagePayment;
            ltvU = max(balancePrime / (pHU * quantity), 0);
            ltvD = max(balancePrime / (pHD * quantity), 0);

            // start the optimization process over a'
            for (int iap = 0; iap < na; iap++){

                // re-writing the budget constraint
                consumption = (1 + r) * agrid[ia] + w * exp(egrid[ie]) - agrid[iap] - mortgagePayment;

                // bequest if house prices go up or down, respectively
                bequestU = (1 + r) * agrid[iap] + (1 - ltvU) * pHU * quantity;
                bequestD = (1 + r) * agrid[iap] + (1 - ltvD) * pHD * quantity;

                /* expectation of the log of bequest (this is important, we are taking expectation over the value of the
                 * bequest, not the value of the expected bequest i.e. E(log(b)) not log(E(b))
                 * (re. Jensen's inequality issues)
                 */
                expected = piU * log(bequestU) + piD * log(bequestD);

                // calculate the utility: flow utility (log/Cobb-Douglas) plus the expected log of bequest
                utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected;

                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV){
                    VV = utility;
                    CC = consumption;
                    QQ = qgrid[iq];
                    AA = agrid[iap];
                }
            }
        }


        double[] output = new double[5];
        output[0] = quantity;
        output[1] = VV;
        output[2] = CC;
        output[3] = AA;
        output[4] = QQ;
        return output;
    }

    double[] value_t(int age, int ia, int iq, int iOLTV, int iLTV, int iFICO, int ie,
                     int na, int nq, int nOLTV,
                     double[] agrid, double[] qgrid, double[] oltvGrid, double[][] creditSurface,
                     double[][] ltvSchedule, double[] egrid, double[][][][][][][] Expected,
                     double alpha, double beta, double w, double r, double pH, int N){

        // mortgage terms
        double mortgageRate = creditSurface[iFICO][iOLTV];              // given the agent's point on the CS, back out the mortgage rate
        double originationLTV = oltvGrid[iOLTV];                        // given the agent's point on the CS, back out the origination ltv
        double quantity = qgrid[iq];                                  // given the agent's point on the quantity grid, back out how much housing they hold
        double existingMortgagePayment =                            // given the current price of housing and how much housing the guy has, calculate the mortgage payment
                mortgageRate * originationLTV * pH * quantity / (1 - pow(1 + mortgageRate, - N));
        double currentLTV = ltvSchedule[iOLTV][iLTV];                  // given how many periods have elapsed ('in'), calculate the current ltv of the mortgage

        double wage = w * exp(egrid[ie]);
        double currentAsset = agrid[ia];

        // in case of refi, sales, default
        double newMortgageRate, newOriginationLTV, newMortgagePayment;
        // new_quantity is needed only in case of buy/sell and default
        double newQuantity;

        // value function terms
        double utility;
        double consumption;
        double defaultPenalty = 1;
        double VV = pow(-10, 5);

        // policy function terms
        double CC = 0;
        double QQ = 0;
        double AA = 0;
        double status = 0;

        // problems for households in their interim periods of life (current, refi, buy/sell)
        double VV_C = pow(-10, 5); //value from staying *C*urrent
        double CC_C = 0;
        double AA_C = 0;
        double VV_R = pow(-10, 5); //value from *R*efinancing
        double CC_R = 0;
        double AA_R = 0;
        double VV_S = pow(-10, 5); // value from *S*elling
        double CC_S = 0;
        double AA_S = 0;
        double QQ_S = 0;
        double VV_D = pow(-10, 5); // value from *D*efault
        double CC_D = 0;
        double AA_D = 0;
        double QQ_D = 0;

        double expected_C, expected_R, expected_S, expected_D;

        /* this is the first main optimization loop: asset choice, other loops for choice variables will be
         * nested within this loop
         * */
        for (int iap = 0; iap < na; iap++) {
            double futureAsset = agrid[iap];
            //#####################################################//
            //--------------STAY CURRENT: START--------------------//
            //#####################################################//
            expected_C = Expected[age + 1][iap][iq][iOLTV][iLTV][iFICO][ie];
            consumption = (1 + r) * currentAsset + wage - existingMortgagePayment - futureAsset;
            utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected_C;
            if (consumption <= 0){utility = pow(-10, 5);}
            if (utility >= VV_C){
                VV_C = utility;
                CC_C = consumption;
                AA_C = futureAsset;
            }
            //#####################################################//
            //--------------STAY CURRENT: END----------------------//
            //#####################################################//

            for (int iOLTVp = 0; iOLTVp < nOLTV; iOLTVp++){
                //#################################################//
                //----------------REFINANCE: START-----------------//
                //#################################################//
                newMortgageRate = creditSurface[iFICO][iOLTVp];
                newOriginationLTV = oltvGrid[iOLTVp];
                newMortgagePayment = newMortgageRate * newOriginationLTV * pH * quantity
                        / (1 - pow(1 + newMortgageRate, -N));
                expected_R = Expected[age + 1][iap][iq][iOLTVp][4][iFICO][ie];
                consumption = (1 + r) * currentAsset + wage + (newOriginationLTV - currentLTV) * pH * quantity
                        - futureAsset - newMortgagePayment;
                utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected_R;
                if (consumption <= 0){utility = pow(-10, 5);}
                if (utility >= VV_R){
                    VV_R = utility;
                    CC_R = consumption;
                    AA_R = futureAsset;
                }
                //##################################################//
                //--------------------REFINANCE: END----------------//
                //##################################################//

                /* Start the nested loop for housing quantity. Buy/sell and and default will be
                 * nested within this loop
                 * */
                for (int iqp = 0; iqp < nq; iqp++) {
                    //##################################################//
                    //---------------BUY/SELL: START--------------------//
                    //##################################################//
                    newQuantity = qgrid[iqp];
                    newMortgagePayment = newMortgageRate * newOriginationLTV * pH * newQuantity
                            / (1 - pow(1 + newMortgageRate, -N));
                    expected_S = Expected[age + 1][iap][iqp][iOLTVp][4][iFICO][ie];
                    consumption = (1 + r) * currentAsset + wage + (1 - currentLTV) * pH * quantity
                            - (1 - newOriginationLTV) * pH * newQuantity - futureAsset - newMortgagePayment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(newQuantity) + beta * expected_S;
                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV_S){
                        VV_S = utility;
                        CC_S = consumption;
                        AA_S = futureAsset;
                        QQ_S = newQuantity;
                    }
                    //#####################################################//
                    //--------------------BUY/SELL: END--------------------//
                    //#####################################################//

                    //#####################################################//
                    //--------------------DEFAULT: START-------------------//
                    //#####################################################//

                    expected_D = expected_S;
                    consumption = (1 + r) * currentAsset + wage - (1 - newOriginationLTV) * pH * newQuantity - futureAsset
                            - newMortgagePayment;
                    utility = (1 - alpha) * log(consumption) + alpha * log(quantity) - defaultPenalty + beta * expected_D;
                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV_D){
                        VV_D = utility;
                        CC_D = consumption;
                        AA_D = futureAsset;
                        QQ_D = newQuantity;
                    }
                    //########################################################//
                    //--------------------DEFAULT: END------------------------//
                    //########################################################//
                }
            }
        }

        VV = max(max(VV_C, VV_S), max(VV_R, VV_D));
        if (VV == VV_C){
            CC = CC_C;
            AA = AA_C;
            QQ = qgrid[iq];
            status = 1;
        }else if (VV == VV_R){
            CC = CC_R;
            AA = AA_R;
            QQ = qgrid[iq];
            status = 2;
        }else if (VV == VV_S){
            CC = CC_S;
            AA = AA_S;
            QQ = QQ_S;
            status = 3;
        }else if (VV == VV_D){
            CC = CC_D;
            AA = AA_D;
            QQ = QQ_D;
            status = 4;
        }


        double[] output = new double[6];
        output[0] = quantity;
        output[1] = VV;
        output[2] = CC;
        output[3] = AA;
        output[4] = QQ;
        output[5] = status;
        return output;
    }

    double[] value_0(int age, int ia, int iFICO, int ie, int na, int nq, int nOLTV,
                     double[] agrid, double[] qgrid, double[] oltvGrid, double[][] creditSurface,
                     double[] egrid, double[][][][][][][] Expected,
                     double alpha, double beta, double w, double pH, int N){

        double currentAsset = agrid[ia];
        double futureAsset;
        double wage = w * egrid[ie];

        double mortgageRate;
        double originationLTV;
        double mortgagePayment;
        double downPayment;
        double quantity = 0;

        double consumption;
        double expected;
        double utility;

        double VV = pow(-10, 5);
        double CC = 0;
        double AA = 0;
        double QQ = 0;
        double status = 5;

        for (int iap = 0; iap < na; iap++){
            futureAsset = agrid[iap];
            for (int iq = 0; iq < nq; iq++){
                quantity = qgrid[iq];
                for (int iOLTV = 0; iOLTV < nOLTV; iOLTV++){
                    mortgageRate = creditSurface[iFICO][iOLTV];
                    originationLTV = oltvGrid[iOLTV];
                    mortgagePayment = mortgageRate * originationLTV * pH * quantity / (1 - pow(1 + mortgageRate, - N));
                    downPayment = (1 - originationLTV) * pH * quantity;
                    expected = Expected[age + 1][iap][iq][iOLTV][4][iFICO][ie];
                    consumption = currentAsset + wage - downPayment - mortgagePayment - futureAsset;
                    utility = (1 - alpha) * log(consumption) + alpha * log(quantity) + beta * expected;
                    if (consumption <= 0){utility = pow(-10, 5);}
                    if (utility >= VV){
                        VV = utility;
                        CC = consumption;
                        AA = futureAsset;
                        QQ = quantity;
                    }
                }
            }
        }

        double[] output = new double[6];
        output[0] = 0;
        output[1] = VV;
        output[2] = CC;
        output[3] = AA;
        output[4] = QQ;
        output[5] = status;
        return output;

    }



}
