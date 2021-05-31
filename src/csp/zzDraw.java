///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package csp;
//
//import java.awt.Color;
//import org.jfree.chart.ChartFactory;
//import org.jfree.chart.ChartPanel;
//import org.jfree.chart.JFreeChart;
//import org.jfree.chart.axis.NumberAxis;
//import org.jfree.chart.plot.CategoryPlot;
//import org.jfree.chart.plot.PlotOrientation;
//import org.jfree.chart.plot.XYPlot;
//import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
//import org.jfree.data.category.CategoryDataset;
//import org.jfree.data.xy.XYDataset;
//import org.jfree.data.xy.XYSeries;
//import org.jfree.data.xy.XYSeriesCollection;
//import org.jfree.ui.ApplicationFrame;
//import org.jfree.ui.RectangleInsets;
//import org.jfree.ui.RefineryUtilities;
//
///**
// *
// * @author Anurag
// */
//
//
//abstract class ScatterGraph extends ApplicationFrame {
///**
//* Creates a new demo.
//*
//* @param title the frame title.
//*/
//
//    public ScatterGraph(String title) {
//        super(title);
//    }
//
//    private JFreeChart createDefaultChartXY(XYDataset dataset) {
//        // create the chart...
//           JFreeChart chart = ChartFactory.createScatterPlot(
//                   "Scatter Demo",
//                    "X", // x axis label
//                    "Y", // y axis label
//                    dataset, // data
//                    PlotOrientation.VERTICAL,
//                    true, // include legend
//                    true, // tooltips
//                    false // urls
//                    );
//        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
//        chart.setBackgroundPaint(Color.GREEN);
//        // get a reference to the plot for further customisation...
//        XYPlot plot = (XYPlot) chart.getPlot();
//        plot.setBackgroundPaint(Color.lightGray);
//        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
//        plot.setDomainGridlinePaint(Color.white);
//        plot.setRangeGridlinePaint(Color.white);
//        
//        XYLineAndShapeRenderer renderer
//        = (XYLineAndShapeRenderer) plot.getRenderer();
//        renderer.setShapesVisible(true);
//        renderer.setShapesFilled(true);
//
//        // change the auto tick unit selection to integer units only...
//        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
//        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
//        // OPTIONAL CUSTOMISATION COMPLETED.
//        return chart;
//    }
//    
//    
//    private JFreeChart createDefaultChartLine(XYDataset dataset) {
//        // create the chart...
//           JFreeChart chart = ChartFactory.createXYLineChart(
//                   "Line Demo",
//                    "X axis label", // x axis label
//                    "Y axis label", // y axis label
//                    dataset, // data
//                    PlotOrientation.VERTICAL,
//                    true, // include legend
//                    true, // tooltips
//                    false // urls
//                    );
//           
//        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
//        chart.setBackgroundPaint(Color.GREEN);
//        // get a reference to the plot for further customisation...
////        CategoryPlot plot = (CategoryPlot)chart.getPlot();
//        XYPlot plot = (XYPlot) chart.getPlot();
//        plot.setBackgroundPaint(Color.lightGray);
//        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
//        plot.setDomainGridlinePaint(Color.white);
//        plot.setRangeGridlinePaint(Color.white);
//        
//        XYLineAndShapeRenderer renderer
//        = (XYLineAndShapeRenderer) plot.getRenderer();
//        renderer.setShapesVisible(true);
//        renderer.setShapesFilled(true);
//
//        // change the auto tick unit selection to integer units only...
//        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
//        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
//        // OPTIONAL CUSTOMISATION COMPLETED.
//        return chart;
//    }    
//    
//    
//    
//    /**
//    * Creates a panel for the demo (used by SuperDemo.java).
//    *
//    * @return A panel.
//    */
////    public JPanel createDemoPanel() {
////        JFreeChart chart = createDefaultChart(createDataset());
////        return new ChartPanel(chart);
////    }
//
//    protected void defaultDrawXY(XYDataset dataset){
//        JFreeChart chart = createDefaultChartXY(dataset);
//        ChartPanel chartPanel = new ChartPanel(chart);
//        //chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
//        setContentPane(chartPanel);
//
//        this.pack();
//        RefineryUtilities.centerFrameOnScreen(this);
//        this.setVisible(true);
//    }
//    
//   protected void defaultDrawLine(XYDataset dataset){
//        JFreeChart chart = createDefaultChartLine(dataset);
//        ChartPanel chartPanel = new ChartPanel(chart);
//        //chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
//        setContentPane(chartPanel);
//
//        this.pack();
//        RefineryUtilities.centerFrameOnScreen(this);
//        this.setVisible(true);
//    }
//}
//
//public class Draw extends ScatterGraph{    
//    Draw(String title){
//        super(title);
//    }
//    
//    private XYDataset updateXYDataset(XYSeriesCollection allseries, XYSeries series){
//
//        XYSeriesCollection dataset = allseries;
//        dataset.addSeries(series);
//
//        return dataset;
//    }    
//    
//    public void drawLine(XYDataset dataset){
//        this.defaultDrawLine(dataset);         
//    }
//    
//    public void drawScatter(XYSeries series, XYSeriesCollection allseries){
//        XYDataset dataset; 
//        dataset = this.updateXYDataset(allseries,series);
//        this.defaultDrawXY(dataset);         
//    }
//}
//  
