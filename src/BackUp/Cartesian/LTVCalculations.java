package BackUp;

import java.util.ArrayList;
import java.util.Collections;

public class LTVCalculations {


    public ArrayList<Double> LTVForecast(double[] input, double originationLTV, double interestRate){

        double originationPrice = 500000; //placeholder
        int introPeriod = 24;
        int term = 360;
        ArrayList<Double> balance = new ArrayList<>();
        ArrayList<Double> HPI = new ArrayList<>();
        ArrayList<Double> HPA = new ArrayList<>();
        ArrayList<Double> HPAPercentChange = new ArrayList<>();
        ArrayList<Double> price = new ArrayList<>();
        ArrayList<Double> LTVForecast = new ArrayList<>();

        double principal = originationPrice * originationLTV;
        double monthlyMortgageRate = interestRate/12;
        double monthlyMortgagePayment = (monthlyMortgageRate * principal)/(1 - Math.pow(1 + monthlyMortgageRate, - term));
        //System.out.println("monthly payment = " + monthlyMortgagePayment);

        //BALANCE SCHEDULE:
        balance.add(0, 0d);
        for (int i = 1; i < term + 1; i++) {
            balance.add(i, (balance.get(i - 1) + monthlyMortgagePayment) / (1 + monthlyMortgageRate));
        }
        Collections.reverse(balance);
        //System.out.println("balance schedule = " + balance);


        for (int i = 0; i < input.length; i++){
            HPI.add(i, input[i]);
        }
        Collections.reverse(HPI);

        ArrayList<Double> weights = new ArrayList<>();
        for (int i = 0; i < input.length; i++) {
            weights.add(i, 0.118926 * Math.pow(2, -(i + 1) / 12d)); //the constant of proportionality is s.t. the weights sum to one, 0.118926
        }

        price.add(0, originationPrice);

        for (int i = 0; i < input.length - 1; i++) {
            HPA.add(i, Math.log(HPI.get(i) / HPI.get(i + 1)));
        }


        for (int j = 0; j < introPeriod; j++) {

            double weightedAverage = 0;
            for (int i = 0; i < input.length - 1; i++) {
                weightedAverage = weightedAverage + (weights.get(i) * HPA.get(i));
            }
            HPA.add(j, weightedAverage);

            double HPAPercent = Math.pow(Math.E, weightedAverage) - 1;
            HPAPercentChange.add(j, HPAPercent);

            price.add(j + 1, price.get(j) * (1 + (HPAPercent)));
        }
        //System.out.println("price = " + price);


        for (int i = 0; i < term - introPeriod ; i++){
            price.add(introPeriod + 1 + i, price.get(introPeriod + i) * (1 + 0.02/12d));
        }
        //System.out.println("price = " + price);

        for (int i = 0; i < price.size(); i++){
            LTVForecast.add(i, balance.get(i)/price.get(i));
        }
        //System.out.println("LTV Forecast = " + LTVForecast);

        return LTVForecast;
    }

}
