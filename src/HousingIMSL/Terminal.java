package HousingIMSL;

import com.imsl.math.NelderMead;
import com.imsl.math.Spline2DInterpolate;

import static java.lang.Math.*;
import static java.lang.Math.max;

public class Terminal {

    private double currentAsset = 0;
    private double quantity = 0;
    private double originationLTV = 0;
    private double fico = 0;
    private double currentLTV = 0;
    private double productivity = 0;
    private Spline2DInterpolate creditSurface;

    public void setCurrentAsset(double currentAsset) {
        this.currentAsset = currentAsset;
    }

    public void setQuantity(double quantity) {
        this.quantity = quantity;
    }

    public void setOriginationLTV(double originationLTV) { this.originationLTV = originationLTV; }

    public void setFico(double fico) { this.fico = fico; }

    public void setCurrentLTV(double currentLTV) {
        this.currentLTV = currentLTV;
    }

    public void setProductivity(double productivity) {
        this.productivity = productivity;
    }

    public void setCreditSurface(Spline2DInterpolate creditSurface) {
        this.creditSurface = creditSurface;
    }

    NelderMead.Function objectiveTerminal = x -> {
        double futureAsset = x[0];

        double mortgageRate = creditSurface.value(fico, originationLTV);

        double alpha = 0.33;
        double beta = 0.97;
        double housePrice = 0.5;
        double housePriceUp = housePrice * exp(0 + 0.06);
        double housePriceDown = housePrice * exp(0 - 0.06);
        double interestRate = 0.03;
        double term = 30;

        double wage = exp(productivity);
        double mortgagePayment = mortgageRate * originationLTV * housePrice * quantity / (1 - pow(1 + mortgageRate, -term));
        if (currentLTV == 0){mortgagePayment = 0;}
        double ltvUp = max(((1 + mortgageRate) * (1.2 * currentLTV) * housePrice * quantity - mortgagePayment) / (housePriceUp * max(quantity, 0.0001)), 0);
        double ltvDown = max(((1 + mortgageRate) * (1.2 * currentLTV) * housePrice * quantity - mortgagePayment) / (housePriceDown * max(quantity, 0.0001)), 0);
        double bequestUp = (1 + interestRate) * futureAsset + (1 - ltvUp) * housePriceUp * quantity;
        double bequestDown = (1 + interestRate) * futureAsset + (1 - ltvDown) * housePriceDown * quantity;
        double expected = 0.5 * log(max(bequestUp, 0.0001)) + 0.5 * log(max(bequestDown, 0.0001));
        double consumption = (1 + interestRate) * currentAsset + wage - mortgagePayment - futureAsset;
        double utility = 0;
        //if (consumption <= 0){utility = -1e5;}
        utility = (1 - alpha) * log(consumption) + alpha * log(max(quantity, 0.0001)) + beta * expected;
        return - utility;
    };
}
