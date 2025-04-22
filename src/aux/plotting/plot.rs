// This module is responsible for plotting data using gnuplot.
extern crate gnuplot;
use gnuplot::{Figure, AxesCommon, Caption, Color, Fix};

pub fn plot(energies: &Vec<f64>, stopping_powers: &Vec<f64>, 
    caption: &str, title: &str){
// --- Plotting using gnuplot ---
let mut fg = Figure::new();
{
    let axes = fg.axes2d();

    // Set logarithmic scales for both axes
    axes.set_x_log(Some(2.0));
    axes.set_y_log(Some(2.0));
    
    // Set axis ranges
    axes.set_x_range(Fix(10.0), Fix(10500.0));
    axes.set_y_range(Fix(1e-30),Fix(1e-27));
    
    // Set titles and labels
    axes.set_title(title, &[]);
    axes.set_x_label("Energía (MeV)", &[]);
    axes.set_y_label("Poder de frenado (MeV/cm)", &[]);
    
    // Plot the data in blue with a label
    axes.lines(energies, stopping_powers, &[Caption(caption), Color("blue")]);
}
// Set terminal to PNG (size 1000x600) and display the plot
fg.set_terminal("pngcairo size 1000,600", &format!("./output/{}.png",title));
fg.show().expect("Unable to show plot");
}