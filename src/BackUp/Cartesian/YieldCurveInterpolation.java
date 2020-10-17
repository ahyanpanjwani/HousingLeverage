package BackUp;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class YieldCurveInterpolation {

    public double[] InterpolatedYieldCurve(double[] yieldCurve){
        double[] interpolatedYC;

        interpolatedYC = new double[361];

        double[] time = {1, 3, 6, 12, 24, 36, 60, 84, 120, 240, 360};

        PolynomialSplineFunction fn = new SplineInterpolator().interpolate(time, yieldCurve);

        for (int i = 1; i < 361; i++){
            interpolatedYC[i] = fn.value(i)/1200;
        }

        //System.out.println(Arrays.toString(interpolatedYC));

        return interpolatedYC;
    }

}

