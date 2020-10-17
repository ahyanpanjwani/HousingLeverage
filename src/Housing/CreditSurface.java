package Housing;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class CreditSurface {

    // grid for ficos
    public double[] ficos(){
        double[] ficos = {0, 1d/6, 2d/6, 3d/6, 4d/6, 5d/6, 1};
        return ficos;
    }

    // grid for oltv (it calls the credit surface csv and makes a grid for the oltvs for fico=500)
    public double[] originationLTVs(String filePathForCS, int nl){
        double[] lgrid = new double[nl];
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(filePathForCS));
            String row;
            int it = 0;
            while ((row = bufferedReader.readLine()) != null & it < nl) {
                String[] data = row.split(",");
                lgrid[it] = Double.parseDouble(data[1]);
                it++;
            }
        }catch (IOException e){e.printStackTrace();}
        return lgrid;
    }

    // matrix for the full credit surface matrix
    public double[][] creditSurface(String filePathForCS, int nf, int nCS){
        List<String[]> rowList2 = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(filePathForCS))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList2.add(lineItems);
            }
            br.close();
        }catch (Exception e){e.printStackTrace();}

        double[][] creditSurface = new double[nf][nCS];

        for (int ifico = 0; ifico < nf; ifico++){
            for (int il = 0; il < nCS; il++){
                creditSurface[ifico][il] = Double.parseDouble(rowList2.get(il + (nCS * ifico))[2]);
            }
        }
        return creditSurface;
    }

}
