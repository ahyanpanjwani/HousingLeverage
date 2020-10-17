package Sparse;

public class ExampleFunction extends sgpp.ScalarFunction {

    private double a = 0;
    private double b = 0;

    public ExampleFunction() {
        super(2);
    }
    public double eval(sgpp.DataVector x) {
        if ((x.get(0) >= 0.0) && (x.get(0) <= 1.0) &&
                (x.get(1) >= 0.0) && (x.get(1) <= 1.0)) {
            // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
            return Math.sin(a * x.get(0)) + Math.sin(b * x.get(1));
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }
    public void clone(sgpp.SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunction_t clone) {
    }

    public void setA(double a){
        this.a = a;
    }
    public void setB(double b){
        this.b = b;
    }
}
