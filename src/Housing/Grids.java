package Housing;

import sgpp.DataVector;
import sgpp.OperationEval;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Grids {

    double[] linspace(int nx, double xmin, double xmax){

        double[] agrid = new double[nx];
        agrid[0] = xmin;
        double astep = (xmax - xmin) / (nx - 1);

        for (int i = 1; i < agrid.length; i++){
            agrid[i] = agrid[i - 1] + astep;
        }

        return agrid;
    }

    void dataWriter(ArrayList<double[]> outputMatrix, String destinationFilePath){

        try(FileWriter fileWriter = new FileWriter(destinationFilePath)){
            String[] header = { "age", "value", "current_asset", "current_quantity",
                                "current_origination_ltv", "fico_score", "current_ltv",
                                "wage", "consumption", "future_asset", "future_quantity",
                                "future_oltv", "status"};
            for (String s : header) {
                fileWriter.append(s);
                fileWriter.append(",");
            }
            fileWriter.append("\n");
            for (double[] row : outputMatrix) {
                for (int j = 0; j < header.length; j++) {
                    fileWriter.append(String.valueOf(row[j]));
                    fileWriter.append(",");
                }
                fileWriter.append("\n");
            }

        }catch (IOException ioe){ioe.printStackTrace();}
    }

    double[] originationLtvGrid(String filePathForCS, int numberOfPoints){
        double[] lgrid = new double[numberOfPoints];
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(filePathForCS));
            String row;
            int it = 0;
            while ((row = bufferedReader.readLine()) != null & it < numberOfPoints) {
                String[] data = row.split(",");
                lgrid[it] = Double.parseDouble(data[1]);
                it++;
            }
        }catch (IOException e){e.printStackTrace();}
        return lgrid;
    }

    double[] expectations(double asset, double quantity, double originationLTV, double ficoScore,
                          double ltvUp, double ltvDown, double[] productivityGrid, double[] transitionVector,
                          DataVector alpha, OperationEval operationEval){

        double[] output = new double[2];
        DataVector pUp = new DataVector(6); DataVector pDown = new DataVector(6);
        pUp.set(0, asset); pUp.set(1, quantity); pUp.set(2, originationLTV); pUp.set(3, ficoScore); pUp.set(4, ltvUp);
        pDown.set(0, asset); pDown.set(1, quantity); pDown.set(2, originationLTV); pDown.set(3, ficoScore); pDown.set(4, ltvDown);

        double expectedUp = 0; double expectedDown = 0;
        for (int iep = 0; iep < productivityGrid.length; iep++){
            double productivity = iep / (double) (productivityGrid.length - 1);
            pUp.set(5, productivity);
            pDown.set(5, productivity);
            expectedUp = expectedUp + transitionVector[iep] * operationEval.eval(alpha, pUp);
            expectedDown = expectedDown + transitionVector[iep] * operationEval.eval(alpha, pDown);
        }
        output[0] = expectedUp;
        output[1] = expectedDown;
        //System.out.println(ltvUp + " , " + ltvDown);
        return output;
    }

}
