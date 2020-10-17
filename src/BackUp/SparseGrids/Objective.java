package BackUp.SparseGrids;

import static java.lang.Math.exp;
import static java.lang.Math.log;

public class Objective extends sgpp.ScalarFunction{

    private double currentAsset = 0;
    private double productivity = 0;
    private double interestRate = 0;
    private double discountRate = 0;

    public Objective(){super(2);}

    public double eval(sgpp.DataVector x){

        sgpp.LoadJSGPPLib.loadJSGPPLib();

        double r = interestRate;
        double beta = discountRate;

        double wage = exp(productivity);
        double futureAsset = x.get(0);
        double bequestU = (1 + 0.04) * futureAsset;
        double bequestD = (1 + 0.02) * futureAsset;
        double expected = 0.5 * log(bequestU) + 0.5 * log(bequestD);
        double consumption = (1 + r) * currentAsset + wage - futureAsset;
        double utility = log(consumption) + beta * expected;
        return - utility;
    }

    public void clone(sgpp.SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunction_t clone) {
    }

    public void setCurrentAsset(double currentAsset){
        this.currentAsset = currentAsset;
    }

    public void setProductivity(double productivity){
        this.productivity = productivity;
    }

    public void setInterestRate(double interestRate){
        this.interestRate = interestRate;
    }

    public void setDiscountRate(double discountRate){
        this.discountRate = discountRate;
    }


}
