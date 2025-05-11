// This module provide the function to calculate the stopping power using the Bethe-Bloch formula
use std::collections::HashMap;
use std::fs::File;
use std::f64::consts::PI;
use std::io::Write;

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

pub fn bethe_bloch_no_corrections(n_points: &u32, energies: &mut Vec<f64>, stopping_powers: &mut Vec<f64>){
    
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
    
    // Open file for writing results
    let mut file = File::create("output/fstopping_no_corrections.txt").expect("Unable to create file");

    // BETHE-BLOCH WITHOUT CORRECTIONS
    println!("Bethe-Bloch without corrections");

    for i in 0..*n_points {
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
}

pub fn bethe_bloch_density_corrections(n_points: &u32, energies: &mut Vec<f64>, stopping_powers: &mut Vec<f64>, variables: &HashMap<String, f64>){
    // Retrieve variables from the HashMap
    let a: f64 = variables.get(&String::from("a")).copied().unwrap();
    let x0: f64 = variables.get(&String::from("x0")).copied().unwrap();
    let x1: f64 = variables.get(&String::from("x1")).copied().unwrap();
    let m_param: f64 = variables.get(&String::from("m_param")).copied().unwrap();
    let c_param: f64 = variables.get(&String::from("c_param")).copied().unwrap();

    // Derived constants
    let proton_mass = (PROTON_MASS_0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / ELECTRON_CHARGE;
    let electron_mass = ELECTRON_MASS_0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    let electron_mass_eV = (ELECTRON_MASS_0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / ELECTRON_CHARGE;

    // General constant for the Bethe–Bloch calculation
    let const_general = (4.0 * PI * ELECTRON_CHARGE.powi(4) * COULOMB_CONST.powi(2)) / (electron_mass * ELECTRON_CHARGE * 1.0e8);

    let mut file = File::create("output/fstopping_density_corrections.txt").expect("Unable to create file");
    
    println!("Bethe-Bloch with Density Corrections");

    for i in 0..*n_points{
        let energy_eV = PROTON_ENERGY_MeV_I * ((i as f64 + 1.0) * 10.0) * 1e6;

        let beta = ((energy_eV * (energy_eV + 2.0 * proton_mass)).sqrt())
                   / (energy_eV + proton_mass);

        let bg = beta * (1.0 / (1.0 - beta * beta).sqrt());

        let x = bg.log10();
        let delta = if x >= x1{
            2.0 * 10.0_f64.log10() * x + c_param   
        } else if x0 <= x {
            2.0 * 10.0_f64.log10() * x + c_param + a * f64::powf(x1 - x,m_param) 
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
}

pub fn bethe_bloch_layer_corrections(n_points: &u32, energies: &mut Vec<f64>, stopping_powers: &mut Vec<f64>){

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


    let mut file = File::create("output/fstopping_layer_corrections.txt").expect("Unable to create file");

    println!("Bethe-Bloch with Layer Correction");

    for i in 0..*n_points{
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
}

pub fn bethe_bloch_all_corrections(n_points: &u32, energies: &mut Vec<f64>, stopping_powers: &mut Vec<f64>, variables: &HashMap<String, f64>) {
    // Retrieve variables from the HashMap
    let a: f64 = variables.get(&String::from("a")).copied().unwrap();
    let x0: f64 = variables.get(&String::from("x0")).copied().unwrap();
    let x1: f64 = variables.get(&String::from("x1")).copied().unwrap();
    let m_param: f64 = variables.get(&String::from("m_param")).copied().unwrap();
    let c_param: f64 = variables.get(&String::from("c_param")).copied().unwrap();

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

    // Open file for writing results
    let mut file = File::create("output/fstopping_all_corrections.txt").expect("Unable to create file");

    // BETHE-BLOCH WITH ALL CORRECTIONS
    println!("Bethe-Bloch with all corrections");

    for i in 0..*n_points{
        let energy_eV = PROTON_ENERGY_MeV_I * ((i as f64 + 1.0) * 10.0) * 1e6;
    
        let beta = ((energy_eV * (energy_eV + 2.0 * proton_mass)).sqrt())
                        / (energy_eV + proton_mass);
            
        let bg = beta * (1.0 / (1.0 - beta * beta).sqrt());

        // delta density correction
        let x = bg.log10();
        let delta = if x >= x1{
            2.0 * 10.0_f64.log10() * x + c_param   
        } else if x0 <= x {
            2.0 * 10.0_f64.log10() * x + c_param + a * f64::powf(x1 - x,m_param) 
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
}