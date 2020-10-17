package Sparse;

import Housing.CreditSurface;
import com.imsl.math.Spline2DInterpolate;
import flanagan.interpolation.BiCubicInterpolation;
import flanagan.interpolation.BiCubicSplineFast;
import flanagan.math.Gradient;
import org.apache.commons.math3.analysis.BivariateFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolator;
import org.apache.commons.math3.analysis.interpolation.BivariateGridInterpolator;

import java.util.Arrays;

public class TestFunction{


    public static void main(String[] args){

        Control control = new Control();

        String filePathForCS = "/home/ahyan/Dropbox/Housing Market and Leverage Cycle/Code/src/Housing/Data/CreditSurface.csv";
        int nFICO = 7, nOLTV = 24;

        Housing.CreditSurface creditSurfaceInterpolation = new CreditSurface();
        double[] ficos = creditSurfaceInterpolation.ficos();
        double[] originationLTVs = creditSurfaceInterpolation.originationLTVs(filePathForCS, nOLTV);
        double[][] surface = creditSurfaceInterpolation.creditSurface(filePathForCS, nFICO, nOLTV);

        //for (int i = 0; i < surface.length; i++){System.out.println(Arrays.toString(surface[i]));}
        System.out.println("=================");

        double fico = 0.5, oltv = 0.5;

        BivariateGridInterpolator bgi = new BicubicInterpolator();
        BivariateFunction func = bgi.interpolate(ficos, originationLTVs, surface);
        double m0 = func.value(fico, oltv);
        System.out.println("m_bgi: " + m0);

        System.out.println("=================");

        BiCubicSplineFast bcs = new BiCubicSplineFast(ficos, originationLTVs, surface);
        double m_bcs = bcs.interpolate(fico, oltv);
        System.out.println("m_bcs: " + m_bcs);

        System.out.println("=================");

        BiCubicInterpolation bci = new BiCubicInterpolation(ficos, originationLTVs, surface, 0);
        double m1 = bci.interpolate(fico, oltv);
        System.out.println("m_bci: " + m1);
        double[][] xDerivatives = bci.getGridDydx1();
        double[][] yDerivatives = bci.getGridDydx2();
        double[][] crossDerivatives = bci.getGridD2ydx1dx2();
        Gradient grad = new Gradient(ficos, originationLTVs, surface);
        double[] grad0 = grad.numDerivAtPoint(fico, oltv);
        System.out.println(Arrays.toString(grad0));
        //for (int i = 0; i < xDerivatives.length; i++){System.out.println(Arrays.toString(xDerivatives[i]));}

        System.out.println("=================");

        Spline2DInterpolate rogue = control.creditSurface(filePathForCS, nFICO, nOLTV);
        double m2 = rogue.value(fico, oltv);
        System.out.println("m_jmsl: " + m2);
    }
}