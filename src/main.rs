mod aux;

use std::collections::HashMap;
use std::env;
use std::io;
use aux::{bethe_bloch, plotting};

extern crate gnuplot;


fn main() {

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

    let mut variables = HashMap::new();
    variables.insert(String::from("a"), a);
    variables.insert(String::from("x0"), x0);
    variables.insert(String::from("x1"), x1);
    variables.insert(String::from("c_param"), c_param);
    variables.insert(String::from("m_param"), m_param);


    // Vectors to store energies (in MeV) and stopping power values (in MeV/cm)
    let mut energies = Vec::with_capacity(1000);
    let mut stopping_powers = Vec::with_capacity(1000);
    let n_points = 1000;


    // BETHE-BLOCH WITHOUT CORRECTIONS

    energies.clear();
    stopping_powers.clear();

    bethe_bloch::bb::bethe_bloch_no_corrections(&n_points, &mut energies, &mut stopping_powers);

    plotting::plot::plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch)", 
    "Poder de Frenado en función de la energía SIN correcciones");


    // BETHE-BLOCH WITH DENSISTY 

    energies.clear();
    stopping_powers.clear();

    bethe_bloch::bb::bethe_bloch_density_corrections(&n_points, &mut energies, &mut stopping_powers, &variables);

    plotting::plot::plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch) Correcion Densidad", 
    "Poder de Frenado en función de la energía con correccion de densidad");


    // BETHE-BLOCH WITH LAYER CORRECTION 

    energies.clear();
    stopping_powers.clear();

    bethe_bloch::bb::bethe_bloch_layer_corrections(&n_points, &mut energies, &mut stopping_powers);

    plotting::plot::plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch) Correcion Capa", 
    "Poder de Frenado en función de la energía con correccion de capa");


    // BETHE-BLOCH WITH ALL CORRECTIONS

    energies.clear();
    stopping_powers.clear();

    bethe_bloch::bb::bethe_bloch_all_corrections(&n_points, &mut energies, &mut stopping_powers, &variables);
    
    plotting::plot::plot(&energies, &stopping_powers, "Protones en Agua (Bethe-Bloch) Correciones Densidad y Capa", 
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
}
