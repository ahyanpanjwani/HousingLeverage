/* THIS FILE IS REDUNDANT!

This is the legacy credit surface with out the stochastic processes for interest rates and HPI.
It interpolates over the yield curve to get the per period interest rate and the HPI is adaptively
extrapolated (weighted avg + 2% for 28 years). I will not be using this anymore
*/

import Jama.Matrix;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class DeterministicCS {

    public Matrix CreditSurface(String MonthYear, double[] HPIYear, double[] yieldCurve) {

        //Input files for probabilities

        String defaultRatesPath = "C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Data\\MonthlyDefaultData\\" +
                "DefaultRates" + "Jan-02" + ".csv";
        String prepayRatesPath = "C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Data\\MonthlyPrepayData\\" +
                "PrepayRates" + "Jan-02" + ".csv";
        String ficoFicoPath =  "C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Data\\FicoFico.csv";


        //Output file

        String outputFileName = "CS" + MonthYear + "CF.csv";

        double[] input = HPIYear;

        LTVCalculations ltvCalculations = new LTVCalculations();
        int numberOfFicoBins = 7;
        int numberOfLTVBins = 12;
        double minLTVBin = 0.3;
        double LTVBinWidth = 0.1;
        int minFicoBin = 500;
        int FicoBinWidth = 50;
        int numberOfInterestGridPoints = 20000-1;
        double stepOfInterestGrid = 0.00001;
        double originationPrice = 500000;
        int term = 360;

        //Interpolated Yield Curve
        YieldCurveInterpolation yieldCurveInterpolation = new YieldCurveInterpolation();
        double[] iCurve = yieldCurveInterpolation.InterpolatedYieldCurve(yieldCurve);

        /////////DEFAULT MATRIX////////////////
        List<String[]> rowList = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(defaultRatesPath))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList.add(lineItems);
            }
            br.close();
        }catch (Exception e){e.printStackTrace();}


        String[][] Default = new String[rowList.size()][];
        for (int i = 0; i < rowList.size(); i++){
            String[] row = rowList.get(i);
            Default[i] = row;
        }


        Matrix DefaultRates = new Matrix(numberOfFicoBins, numberOfLTVBins);
        for (int j = 0; j < DefaultRates.getRowDimension(); j++) {
            for (int i = 0; i < DefaultRates.getColumnDimension(); i++) {
                DefaultRates.set(j, i, Double.parseDouble(Default[j+1][i + 1]));
            }
        }

        //Loop through default matrix to print:
        //double[][] defaultArrays = DefaultRates.getArray();
        //for (int i = 0; i < DefaultRates.getRowDimension(); i++){ System.out.println(Arrays.toString(defaultArrays[i]));}
        //System.out.println("CUT");

        /////////PREPAY MATRIX////////////////
        List<String[]> rowList2 = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(prepayRatesPath))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList2.add(lineItems);
            }
            br.close();
        }catch (Exception e){}


        String[][] Prepay = new String[rowList2.size()][];
        for (int i = 0; i < rowList2.size(); i++){
            String[] row = rowList2.get(i);
            Prepay[i] = row;
        }



        Matrix PrepayRates = new Matrix(numberOfFicoBins, numberOfLTVBins);
        for (int j = 0; j < PrepayRates.getRowDimension(); j++) {
            for (int i = 0; i < PrepayRates.getColumnDimension(); i++) {
                PrepayRates.set(j, i, Double.parseDouble(Prepay[j+1][i + 1]));
            }
        }

        //Loop through prepay matrix to print:
        //double[][] prepayArrays = PrepayRates.getArray();
        //for (int i = 0; i < PrepayRates.getRowDimension(); i++){ System.out.println(Arrays.toString(prepayArrays[i]));}
        //System.out.println("CUT");

        /////////FICOFICO MATRIX////////////////
        List<String[]> rowList3 = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(ficoFicoPath))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList3.add(lineItems);
            }
            br.close();
        }catch (Exception e){}


        String[][] FICOFICO = new String[rowList3.size()][];
        for (int i = 0; i < rowList3.size(); i++){
            String[] row = rowList3.get(i);
            FICOFICO[i] = row;
        }


        Matrix FicoFicoRates = new Matrix(numberOfFicoBins, numberOfFicoBins);
        for (int j = 0; j < FicoFicoRates.getRowDimension(); j++) {
            for (int i = 0; i < FicoFicoRates.getColumnDimension(); i++) {
                FicoFicoRates.set(j, i, Double.parseDouble(FICOFICO[j+1][i + 1]));
            }
        }

        //Loop through fico fico matrix to print:
        //double[][] FicoFicoArrays = FicoFicoRates.getArray();
        //for (int i = 0; i < FicoFicoRates.getRowDimension(); i++){ System.out.println(Arrays.toString(FicoFicoArrays[i]));}

        //LTV Array
        ArrayList<Double> LTVArray = new ArrayList<>();
        LTVArray.add(0, minLTVBin);
        for (int i = 1; i < numberOfLTVBins; i++){LTVArray.add(i, LTVArray.get(i - 1) + LTVBinWidth);}
        //System.out.println("LTV Array = " + LTVArray);

        //Fico Array
        ArrayList<Integer> FicoArray = new ArrayList<>();
        FicoArray.add(0, minFicoBin);
        for (int i = 1; i < numberOfFicoBins; i++){FicoArray.add(i, FicoArray.get(i - 1) + FicoBinWidth);}
        //System.out.println("FICO Array = " + FicoArray);

        //Interest Rate Grids
        ArrayList<Double> interestRateGrid = new ArrayList<>();
        interestRateGrid.add(0, 0.00001);
        for (int i = 1; i < numberOfInterestGridPoints + 1; i++){interestRateGrid.add(i, interestRateGrid.get(i - 1) + stepOfInterestGrid);}
        //System.out.println("initial grid point = " + interestRateGrid.get(0) + " terminal grid point = " + interestRateGrid.get(numberOfInterestGridPoints));

        Matrix differenceMatrix = new Matrix(numberOfFicoBins, interestRateGrid.size());
        Matrix creditSurface = new Matrix(numberOfLTVBins - 4, numberOfFicoBins);

        for (int l = 0; l < LTVArray.size() - 4; l++) {

            double originationLTV = LTVArray.get(l);
            //System.out.println("LTV = " + originationLTV);
            double principal = originationLTV * originationPrice;


            for (int r = 0; r < interestRateGrid.size(); r++) {
                double interestRate = interestRateGrid.get(r);
                double monthlyInterestRate = interestRate / 12;
                double monthlyMortgagePayment = monthlyInterestRate * principal / (1 - Math.pow(1 + monthlyInterestRate, -term));
                ArrayList<Double> balance = new ArrayList<>();
                balance.add(0, 0d); //terminal balance equals zero
                //recursively calculate the entire balance schedule
                for (int i = 1; i < term + 1; i++) {
                    balance.add(i, (balance.get(i - 1) + monthlyMortgagePayment) / (1 + monthlyInterestRate));
                }

                Collections.reverse(balance);

                ArrayList<Double> ltvForecast = ltvCalculations.LTVForecast(input, originationLTV, interestRate);

                ////////VALUE CALCULATION//////////////
                ArrayList<Double> valueFuture = new ArrayList<>();
                ArrayList<Double> valueCurrent = new ArrayList<>();


                for (int i = 0; i < numberOfFicoBins; i++) {
                    double recoveryFraction = 1;
                    valueCurrent.add(i, ((monthlyMortgagePayment * (1 - DefaultRates.get(i, 0) - PrepayRates.get(i, 0))) / (1 + iCurve[term])) + DefaultRates.get(i, l) * recoveryFraction * balance.get(0) + PrepayRates.get(i, l) * balance.get(0));
                }
                //System.out.println(valueCurrent);

                int t = term;
                while (t > 0) {
                    t--;
                    double recoveryFraction = Math.min(1, 0.5/ltvForecast.get(t));
                    int k = 0;
                    if (ltvForecast.get(t) <= 0.3){k = 0;}
                    else if (ltvForecast.get(t) > 0.3 && ltvForecast.get(t) <= 0.4){k = 1;}
                    else if (ltvForecast.get(t) > 0.4 && ltvForecast.get(t) <= 0.5){k = 2;}
                    else if (ltvForecast.get(t) > 0.5 && ltvForecast.get(t) <= 0.6){k = 3;}
                    else if (ltvForecast.get(t) > 0.6 && ltvForecast.get(t) <= 0.7){k = 4;}
                    else if (ltvForecast.get(t) > 0.7 && ltvForecast.get(t) <= 0.8){k = 5;}
                    else if (ltvForecast.get(t) > 0.8 && ltvForecast.get(t) <= 0.9){k = 6;}
                    else if (ltvForecast.get(t) > 0.9 && ltvForecast.get(t) <= 1.0){k = 7;}
                    else if (ltvForecast.get(t) > 1.0 && ltvForecast.get(t) <= 1.1){k = 8;}
                    else if (ltvForecast.get(t) > 1.1 && ltvForecast.get(t) <= 1.2){k = 9;}
                    else if (ltvForecast.get(t) > 1.2 && ltvForecast.get(t) <= 1.3){k = 10;}
                    else if (ltvForecast.get(t) > 1.3){k = 11;}


                    valueFuture = valueCurrent;
                    for (int i = 0; i < numberOfFicoBins; i++) {
                        valueCurrent.set(i, (((FicoFicoRates.get(i, 0) * (monthlyMortgagePayment + valueFuture.get(0)) + FicoFicoRates.get(i, 1) * (monthlyMortgagePayment + valueFuture.get(1))
                                + FicoFicoRates.get(i, 2) * (monthlyMortgagePayment + valueFuture.get(2)) + FicoFicoRates.get(i, 3) * (monthlyMortgagePayment + valueFuture.get(3))
                                + FicoFicoRates.get(i, 4) * (monthlyMortgagePayment + valueFuture.get(4)) + FicoFicoRates.get(i, 5) * (monthlyMortgagePayment + valueFuture.get(5))
                                + FicoFicoRates.get(i, 6) * (monthlyMortgagePayment + valueFuture.get(6))) * (1 - DefaultRates.get(i, k) - PrepayRates.get(i, k))) / (1 + iCurve[t]))
                                + DefaultRates.get(i, k) * recoveryFraction * balance.get(t) + PrepayRates.get(i, k) * balance.get(t));
                    }
                    //System.out.println(valueCurrent.get(0));
                    //System.out.println(t);
                }
                //System.out.println("Initial Loan Value = " + valueCurrent);
                //System.out.println("principal = " + principal);
                for (int p = 0; p < numberOfFicoBins; p++) {
                    differenceMatrix.set(p, r, Math.abs(valueCurrent.get(p) - principal));
                }
            }

            double couponRate = 0;
            //double[][] differenceArray = differenceMatrix.getArray();
            for (int j = 0; j < numberOfFicoBins; j++) {
                double min = differenceMatrix.get(j, 0);
                int index = 0;
                for (int i = 0; i < differenceMatrix.getColumnDimension(); i++) {
                    if (differenceMatrix.get(j, i) < min) {
                        min = differenceMatrix.get(j, i);
                        index = i;
                    }
                }
                //System.out.println("min difference = " + min + " index = " + index);
                couponRate = interestRateGrid.get(index);
                creditSurface.set(l, j, couponRate * 100);
                //System.out.println("Coupon Rate = " + creditSurface.get(l,j));
            }
        }


//Print Credit Surface Matrix:
        Matrix creditSurfaceTranspose = creditSurface.transpose();
        double[][] creditSurfaceArray = creditSurfaceTranspose.getArray();
        //for (int i = 0; i < DefaultRates.getRowDimension(); i++){System.out.println(Arrays.toString(creditSurfaceArray[i]));}

        ///////////////INTERPOLATION////////////////
        int interpolationPoints = 36;
        double interpolationStep = 0.02;

        Matrix interpolatedCreditSurface = new Matrix(numberOfFicoBins, interpolationPoints);
        ArrayList<Double> interpolationLTV = new ArrayList<>();
        for (int i = 0; i < interpolationPoints; i++){
            interpolationLTV.add(i, minLTVBin + (interpolationStep * i));
        }

        double[] LTV = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

        UnivariateInterpolator univariateInterpolator = new SplineInterpolator();

        for (int j = 0; j < numberOfFicoBins; j++) {
            UnivariateFunction univariateFunction = univariateInterpolator.interpolate(LTV, creditSurfaceArray[j]);
            double interpolatedInterestRate = 0;
            for (int i = 0; i < interpolationPoints; i++) {
                interpolatedInterestRate = Math.max(0, univariateFunction.value(interpolationLTV.get(i)));
                interpolatedCreditSurface.set(j, i, interpolatedInterestRate);
            }
        }


        try(FileWriter fileWriter = new FileWriter("C:\\Users\\ahyan\\Dropbox\\CreditSurfaceTheory\\Output\\CSCounterfactualTwoPercent\\" + outputFileName)){
            for (int i = 0; i < interpolatedCreditSurface.getRowDimension(); i++){
                for (int j = 0; j < interpolatedCreditSurface.getColumnDimension(); j++){
                    fileWriter.append(String.valueOf(FicoArray.get(i)));
                    fileWriter.append(",");
                    fileWriter.append(String.valueOf(interpolationLTV.get(j)));
                    fileWriter.append(",");
                    fileWriter.append(String.valueOf(interpolatedCreditSurface.get(i, j)));
                    //if (j < creditSurface.getRowDimension() - 1){fileWriter.append(",");}
                    fileWriter.append("\n");
                }
            }
            fileWriter.close();
        }catch (IOException ioe){ioe.printStackTrace();}




        return interpolatedCreditSurface;


    }

}