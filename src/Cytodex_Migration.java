
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.filter.BackgroundSubtracter;
import ij.plugin.filter.DifferenceOfGaussians;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.frame.RoiManager;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import java.awt.Color;
import java.awt.Polygon;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class Cytodex_Migration implements PlugIn {

    
     private final boolean canceled =false;
     private boolean isStack = false;
     private Roi cytodexRoi;
     
    // clear outside selection with true background
    private void clearOutside(ImagePlus img, Roi cropRoi) {
        ImageProcessor ip = img.getProcessor();
        ip.setColor(Color.BLACK);
        img.setActivated();
        for(int z = 1; z <= img.getNSlices(); z++) {
            img.setSlice(z);
            ip.fillOutside(cropRoi);
        }
        img.updateAndDraw();
        img.deleteRoi();
    }
    
     // Substract background
    private void backgroundSubstract(ImagePlus img) {
        BackgroundSubtracter imgSubstract = new BackgroundSubtracter();
        if (isStack) {
            for (int s = 1;s <= img.getNSlices(); s++) { 
                img.setSlice(s);
                imgSubstract.rollingBallBackground(img.getProcessor(), 50, false, false, false, false, false);
            }
            img.updateAndRepaintWindow();
        }
        else imgSubstract.rollingBallBackground(img.getProcessor(), 50, false, false, false, false, false);
    }
    
    /**
     * Find nucleus and return the distance to the cytodex
     */
    private double[] findNucleus(ImagePlus imgCrop, String dir , String fileName, int cytodex) {
	int xCoor, yCoor;

        imgCrop.deleteRoi();
        ImagePlus imgNucleus = new Duplicator().run(imgCrop, 1, 1);
        ImagePlus colorNucleus = imgNucleus.duplicate();
        ImageConverter imgConv = new ImageConverter(colorNucleus);
        imgConv.convertToRGB();

    // run difference of Gaussians
        double sigma1 = 4, sigma2 = 1;
        DifferenceOfGaussians.run(imgNucleus.getProcessor(), sigma1, sigma2);
        IJ.run("Contrast Enhancer");
        imgNucleus.getProcessor().setColor(Color.BLACK);
        imgNucleus.setRoi(cytodexRoi);
        imgNucleus.getProcessor().fill(imgNucleus.getRoi());
        imgNucleus.deleteRoi();
// find maxima
        ImageProcessor ipNucleus = imgNucleus.getProcessor();
        MaximumFinder findMax = new MaximumFinder();
        Polygon  nucleusPoly = findMax.getMaxima(ipNucleus, 2,false);        
        double[] nucleusDist = new double[nucleusPoly.npoints];
        colorNucleus.setColor(Color.BLUE);
        for (int i = 0; i < nucleusPoly.npoints; i++) {
            xCoor = (int)nucleusPoly.xpoints[i];
            yCoor = (int)nucleusPoly.ypoints[i];
            nucleusDist[i] = getClosestDistance(xCoor, yCoor, cytodexRoi.getPolygon())*imgCrop.getCalibration().pixelWidth;
            OvalRoi ptNucleus = new OvalRoi(xCoor, yCoor,4,4);
            colorNucleus.getProcessor().draw(ptNucleus);
            colorNucleus.updateAndDraw();                  
        }

        FileSaver imgSave = new FileSaver(colorNucleus);
        imgSave.saveAsTiff(dir+fileName+"_Crop_"+cytodex+"_nucleus.tif");
        colorNucleus.close();
        imgNucleus.close();
        imgNucleus.flush();
        imgCrop.changes = false;
        imgCrop.close();
        imgCrop.flush();
        return nucleusDist;
    }
    
    /**
     * Find nucleus distance to cytodex
     * @param x nucleus positions
     * @param y nucleus positions
     * @param points roi
     * @return 
     */
    double getClosestDistance(int x, int y, Polygon points) {
        double distance = Double.MAX_VALUE;
        for (int i = 0; i < points.npoints; i++) {
                double dx = points.xpoints[i] - x;
                double dy = points.ypoints[i] - y;
                double distance2 = Math.sqrt(dx*dx+dy*dy);
                if (distance2 < distance) {
                        distance = distance2;
                }
        }
        return distance;
    }
             
    /**
     * @param args the command line arguments
     */
    public void run(String arg) {
        if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
        }
        
        // Ask for images folder
        
        String imageDir = IJ.getDirectory("Choose Directory Containing TIF Files...");
        if (imageDir == null) return;
        File inDir = new File(imageDir);
        String [] imageFile = inDir.list();
        if (imageFile == null) return;
        
        // create directory to store images
        String imgOutDir = imageDir+"Images/";
        File imgTmpDir = new File(imgOutDir);
        if (!imgTmpDir.isDirectory())
            imgTmpDir.mkdir();
        
        // write headers for global results
        try {
            FileWriter fwAnalyze;
         
            fwAnalyze = new FileWriter(imageDir + "Analyze_distance_results.xls",false);
         
            BufferedWriter outputAnalyze = new BufferedWriter(fwAnalyze);
            outputAnalyze.write("Image\t#Cytodex\tDistance(µm)\n");
            outputAnalyze.flush();
       
        Duplicator imgDup = new Duplicator();
        
        // Open all images tif in inDir folder  
        for (int i = 0; i < imageFile.length; i++) {
            if (imageFile[i].endsWith(".tif")) {
                String imagePath = imageDir + imageFile[i];
                Opener imgOpener = new Opener();
                ImagePlus imgOrg = imgOpener.openImage(imagePath);
                if (imgOrg.getNSlices() > 1) {
                        isStack = true;
                }
                String fileNameWithOutExt = imageFile[i].substring(0, imageFile[i].length() - 4);
        
        //convert to 8 bytes
                if (isStack && imgOrg.getBytesPerPixel() == 2) {
                    new StackConverter(imgOrg).convertToGray8();
                }
                else if (imgOrg.getBytesPerPixel() == 2) {
                    new ImageConverter(imgOrg).convertToGray8();
                }
           
                imgOrg.show();
        // auto contrast                    
                for ( int s = 1; s <= imgOrg.getNSlices(); s++) {
                    imgOrg.setSlice(s);
                    IJ.run("Contrast Enhancer");   
                }
                if (RoiManager.getInstance() != null) RoiManager.getInstance().close();
                RoiManager rm = new RoiManager();
                new WaitForUserDialog("Select part(s) of image to analyze\nPress t to add selection to the ROI manager.").show();
                
        // find roi 
                for (int r = 0; r < rm.getCount(); r++) {
                    WindowManager.setTempCurrentImage(imgOrg);
                    rm.select(r);
                    int [] roiType = new int[rm.getCount()];
                    roiType[r] = imgOrg.getRoi().getType();
                    IJ.run("Duplicate...", "title=Crop duplicate range=1-" + imgOrg.getNSlices());
                    ImagePlus imgCrop = WindowManager.getCurrentImage();
                    imgCrop.show();
                    
        // substract background                    
                    backgroundSubstract(imgCrop);
                    
        // fill outside roi with background except if ROI is a rectangle
                if (roiType[r] != 0) clearOutside(imgCrop, imgCrop.getRoi()); 
                        
        // save croped image
                FileSaver imgCrop_save = new FileSaver(imgCrop);
                if (isStack) 
                    imgCrop_save.saveAsTiffStack(imgOutDir+fileNameWithOutExt+"_Crop"+r+".tif");
                else 
                    imgCrop_save.saveAsTiff(imgOutDir+fileNameWithOutExt+"_Crop"+r+".tif"); 
                
                if (isStack) 
                    imgCrop.setSlice(2);
                //IJ.run(imgOrg, "Enhance Contrast", "saturated=0.35 normalize"); 
                IJ.run("Contrast Enhancer");
        // ask for outline cytodex
                new WaitForUserDialog("Outline the cytodex").show(); 			
                cytodexRoi = imgCrop.getRoi();

        // find nucleus positions and compute distance to cytodex
                
                double [] dist = findNucleus(imgCrop, imgOutDir, fileNameWithOutExt, r);
                
        // print distance from nucleus to cytodex
            
                for (int d = 0; d < dist.length; d++) {
                // write data
                    outputAnalyze.write(fileNameWithOutExt + "\t" + r + "\t" + dist[d] + "\n");
                    outputAnalyze.flush();
                }   
            }             
            imgOrg.close();
            imgOrg.flush();
            }   
        }
        outputAnalyze.close();
        IJ.showStatus("End of process");
    } catch (IOException ex) {
            Logger.getLogger(Cytodex_Migration.class.getName()).log(Level.SEVERE, null, ex);
    }    
    }
    
}
