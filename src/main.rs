use std::collections::btree_map;
use std::env;
use std::io::{self, Write};
use std::fs::File;
use std::f64::consts::PI;
use std::iter::Enumerate;

extern crate gnuplot;
use gnuplot::{Figure, AxesCommon, Caption, Color, Fix};


fn main() {

    // Physical constants (SI units and energy in eV unless noted)
    const ELECTRON_CHARGE: f64= 1.602176634e-19;
    const ELECTRON_MASS_0: f64 = 9.10938356e-31;
    const SPEED_OF_LIGHT: f64 = 299792458.0;
    const WATER_ATOMIC_NUMBER: f64 = 9.0;
    const PROTON_ENERGY_MeV_I: f64 = 1.0;
    const PROTON_MASS_0: f64 = 1.6726219e-27;
    const WATER_EXCITATION_ENERGY: f64 = 74.6;
    const ELECTRON_PER_VOLUME_H20: f64 = 3.3429;
    const COULOMB_CONST: f64 = 8.99e9;
    const Z_PROTON: f64 = 1.0;

    // Default delta correction parameters
    let a_default = 0.09116;
    let x0_default = 0.24;
    let x1_default = 2.8004;
    let c_default = 3.5017;
    let m_default = 3.4773;

    // Variables for delta correction parameters (initialize with defaults)
    let mut a = a_default;
    let mut x0 = x0_default;
    let mut x1 = x1_default;
    let mut c_param = c_default;
    let mut m_param = m_default;

    // Process command-line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() >= 6 {
        a = args[1].parse().unwrap_or(a_default);
        x0 = args[2].parse().unwrap_or(x0_default);
        x1 = args[3].parse().unwrap_or(x1_default);
        c_param = args[4].parse().unwrap_or(c_default);
        m_param = args[5].parse().unwrap_or(m_default);
    } else {
        println!("Delta correction parameters were not fully provided on the command line.");
        println!("Would you like to input them via standard input? (y/n): ");
        let mut answer = String::new();
        io::stdin().read_line(&mut answer).expect("Failed to read line");
        if answer.trim().eq_ignore_ascii_case("y") {
            a = prompt("Enter value for a", a_default);
            x0 = prompt("Enter value for x0", x0_default);
            x1 = prompt("Enter value for x1", x1_default);
            c_param = prompt("Enter value for C", c_default);
            m_param = prompt("Enter value for m", m_default);
        } else {
            println!("Using default delta correction parameters.");
        }
    }

    // Derived constants
    let proton_mass = (PROTON_MASS_0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / ELECTRON_CHARGE;
    let electron_mass = ELECTRON_MASS_0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    let electron_mass_eV = (ELECTRON_MASS_0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / ELECTRON_CHARGE;

    // General constant for the Bethe–Bloch calculation
    let const_general = (4.0 * PI * ELECTRON_CHARGE.powi(4) * COULOMB_CONST.powi(2)) / (electron_mass * ELECTRON_CHARGE * 1.0e8);

    // Ionization constant I (not used further in the calculation)
    let I = if WATER_ATOMIC_NUMBER < 13.0 {
        (12.0 * WATER_ATOMIC_NUMBER + 7.0) / 1e6
    } else {
        (9.76 * Z_PROTON + 58.8 * WATER_ATOMIC_NUMBER.powf(-0.19)) / 1e6
    };

    // Vectors to store energies (in MeV) and stopping power values (in MeV/cm)
    let mut energies = Vec::with_capacity(1000);
    let mut stopping_powers = Vec::with_capacity(1000);

    // // Open file for writing results
    let mut file = File::create("fstopping_no_corrections.txt").expect("Unable to create file");

    // BETHE-BLOCH WITHOUT CORRECTIONS
    println!("Bethe-Bloch without corrections");

    let n_points = 1000;
    for i in 0..n_points {
        // Calculate proton energy in eV
        let energy_eV = PROTON_ENERGY_MeV_I * ((i as f64 + 1.0) * 10.0) * 1e6;

        // Calculate beta (v/c)
        // beta = sqrt(E*(E + 2*m)) / (E + m)
        let beta = ((energy_eV * (energy_eV + 2.0 * proton_mass)).sqrt())
                   / (energy_eV + proton_mass);

        // Bethe-Bloch formula for the stopping power (dE/dx)
        let de_dx = (const_general * Z_PROTON.powi(2) * ELECTRON_PER_VOLUME_H20 / (beta * beta))
                    * ((2.0 * electron_mass_eV * beta * beta / WATER_EXCITATION_ENERGY).ln()
                       - (1.0 - beta * beta).ln()
                       - beta * beta);

        let energy_MeV = energy_eV / 1e6;

        energies.push(energy_MeV);
        stopping_powers.push(de_dx);

        // Write to file
        writeln!(file, "{:.1}\t{:e}", energy_MeV, de_dx).expect("Unable to write data");
        
        println!("{:.1} MeV (dE/dx): {} MeV/cm", energy_MeV, de_dx);
    }

    plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch)", 
    "Poder de Frenado en función de la energía SIN correcciones");

    // BETHE-BLOCH WITH DENSISTY 

    energies.clear();
    stopping_powers.clear();

    file = File::create("fstopping_density_corrections.txt").expect("Unable to create file");
    
    println!("Bethe-Bloch with Density Corrections");
    
    for i in 0..n_points{
        let energy_eV = PROTON_ENERGY_MeV_I * ((i as f64 + 1.0) * 10.0) * 1e6;

        let beta = ((energy_eV * (energy_eV + 2.0 * proton_mass)).sqrt())
                   / (energy_eV + proton_mass);

        let bg = beta * (1.0 / (1.0 - beta * beta).sqrt());

        let x = bg.log10();
        let delta = if x >= x1{
            2.0 * 10.0_f64.log10() * x - c_param // TODO: Revisar el parametro c_param, debería calcularse y no restarse ?  
        } else if x0 <= x {
            2.0 * 10.0_f64.log10() * x - c_param + a * f64::powf(x1 - x,m_param) 
        } else{
            0.0_f64
        };

        let de_dx = ((const_general * Z_PROTON.powi(2) * ELECTRON_PER_VOLUME_H20) / (beta.powi(2)))
                * ((2.0 * electron_mass_eV * beta.powi(2) / WATER_EXCITATION_ENERGY).ln()
                - (1.0 - beta.powi(2)).ln() - beta.powi(2) - delta);

        let energy_MeV = energy_eV / 1e6;

        energies.push(energy_MeV);
        stopping_powers.push(de_dx);

        // Write to file
        writeln!(file, "{:.1}\t{:e}", energy_MeV, de_dx).expect("Unable to write data");
        println!("{:.1} MeV (dE/dx): {} MeV/cm", energy_MeV, de_dx);
    }

    plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch) Correcion Densidad", 
    "Poder de Frenado en función de la energía con correccion de densidad");

    // BETHE-BLOCH WITH LAYER CORRECTION 

    energies.clear();
    stopping_powers.clear();

    file = File::create("fstopping_layer_corrections.txt").expect("Unable to create file");

    println!("Bethe-Bloch with Layer Correction");

    for i in 0..n_points{
        let energy_eV = PROTON_ENERGY_MeV_I * ((i as f64 + 1.0) * 10.0) * 1e6;

        let beta = ((energy_eV * (energy_eV + 2.0 * proton_mass)).sqrt())
                    / (energy_eV + proton_mass);
        
        let bg = beta * (1.0 / (1.0 - beta * beta).sqrt());

        // shell correction
        let sc = (0.422377*bg.powi(-2) + 0.0304043*bg.powi(-4) - 0.00038106*bg.powi(-6))*(10.0_f64.powi(-6))*((I.powi(2)*10.0_f64.powi(-6))) 
                    + (3.850190*bg.powi(-2)-0.1667989*bg.powi(-4) + 0.00157955*bg.powi(-6))*(10.0_f64.powi(-9))*((I.powi(3)*10.0_f64.powi(-6)));

        let de_dx = ((const_general * Z_PROTON.powi(2) * ELECTRON_PER_VOLUME_H20) / (beta.powi(2)))
        * ((2.0 * electron_mass_eV * beta.powi(2) / WATER_EXCITATION_ENERGY).ln()
        - (1.0 - beta.powi(2)).ln() - beta.powi(2) - 2.0*(sc / WATER_ATOMIC_NUMBER));

        let energy_MeV = energy_eV / 1e6;

        energies.push(energy_MeV);
        stopping_powers.push(de_dx);

            
        writeln!(file, "{:.1}\t{:e}", energy_MeV, de_dx).expect("Unable to write data");
        println!("{:.1} MeV (dE/dx): {} MeV/cm", energy_MeV, de_dx);
    }

    plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch) Correcion Capa", 
    "Poder de Frenado en función de la energía con correccion de capa");

        // BETHE-BLOCH WITH ALL CORRECTIONS

        energies.clear();
        stopping_powers.clear();
    
        file = File::create("fstopping_all_corrections.txt").expect("Unable to create file");
    
        println!("Bethe-Bloch with All Correction");
    
        for i in 0..n_points{
            let energy_eV = PROTON_ENERGY_MeV_I * ((i as f64 + 1.0) * 10.0) * 1e6;
    
            let beta = ((energy_eV * (energy_eV + 2.0 * proton_mass)).sqrt())
                        / (energy_eV + proton_mass);
            
            let bg = beta * (1.0 / (1.0 - beta * beta).sqrt());

            // delta density correction
            let x = bg.log10();
            let delta = if x >= x1{
                2.0 * 10.0_f64.log10() * x - c_param // TODO: Revisar el parametro c_param, debería calcularse y no restarse ?  
            } else if x0 <= x {
                2.0 * 10.0_f64.log10() * x - c_param + a * f64::powf(x1 - x,m_param) 
            } else{
                0.0_f64
            };
    
            // shell correction
            let sc = (0.422377*bg.powi(-2) + 0.0304043*bg.powi(-4) - 0.00038106*bg.powi(-6))*(10.0_f64.powi(-6))*((I.powi(2)*10.0_f64.powi(-6))) 
                        + (3.850190*bg.powi(-2)-0.1667989*bg.powi(-4) + 0.00157955*bg.powi(-6))*(10.0_f64.powi(-9))*((I.powi(3)*10.0_f64.powi(-6)));
    
            let de_dx = ((const_general * Z_PROTON.powi(2) * ELECTRON_PER_VOLUME_H20) / (beta.powi(2)))
            * ((2.0 * electron_mass_eV * beta.powi(2) / WATER_EXCITATION_ENERGY).ln()
            - (1.0 - beta.powi(2)).ln() - beta.powi(2) - delta - 2.0*(sc / WATER_ATOMIC_NUMBER));
    
            let energy_MeV = energy_eV / 1e6;
    
            energies.push(energy_MeV);
            stopping_powers.push(de_dx);
    
                
            writeln!(file, "{:.1}\t{:e}", energy_MeV, de_dx).expect("Unable to write data");
            println!("{:.1} MeV (dE/dx): {} MeV/cm", energy_MeV, de_dx);
        }
    
        plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch) Correciones Densidad y Capa", 
        "Poder de Frenado en función de la energía con correcciones de densidad y capa");

}

// Helper function to prompt the user for a value with a default.
fn prompt(message: &str, default: f64) -> f64 {
    println!("{} (default {}): ", message, default);
    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");

    let trimmed = input.trim();          
    if trimmed.is_empty() {
       default
    } else {
        trimmed.parse().unwrap_or(default)
    }

    // My implementation 
    // match input.trim().parse(){
    //     Ok(num) => num,
    //     Error(_) => default
    // }
}

fn plot(energies: &Vec<f64>, stopping_powers: &Vec<f64>, 
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
    fg.set_terminal("pngcairo size 1000,600", &format!("{}.png",title));
    fg.show().expect("Unable to show plot");
}
