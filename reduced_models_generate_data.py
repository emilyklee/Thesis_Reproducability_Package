import numpy
import copy
import ast
import time
import sys
from scipy.integrate import solve_ivp
from PythonBFM.BFM_reduced_rate_eqns import bfm_reduced_rate_eqns

if __name__ == '__main__':
    
    # Names of species in the system
    species_names = ['O2o', 'N1p', 'N3n', 'N4n', 'O4n', 'N5s', 'N6r', 'B1c', 'B1n', 'B1p', 
                     'P1c', 'P1n', 'P1p', 'P1l', 'P1s', 'P2c', 'P2n', 'P2p', 'P2l',
                     'P3c', 'P3n', 'P3p', 'P3l', 'P4c', 'P4n', 'P4p', 'P4l',
                     'Z3c', 'Z3n', 'Z3p', 'Z4c', 'Z4n', 'Z4p', 'Z5c', 'Z5n', 'Z5p',
                     'Z6c', 'Z6n', 'Z6p', 'R1c', 'R1n', 'R1p', 'R2c', 'R3c', 'R6c', 
                     'R6n', 'R6p', 'R6s', 'O3c', 'O3h']
    
    # Initial concentrations
    c0 = [300.0,                    # O2o
          1.0,                      # N1p
          5.0,                      # N3n
          1.0,                      # N4n
          200.0,                    # O4n
          8.0,                      # N5s
          1.0,                      # N6r
          1.0,                      # B1c
          1.67e-2,                  # B1n
          1.85e-3,                  # B1p
          1.0,                      # P1c
          1.26e-2,                  # P1n
          7.86e-4,                  # P1p
          2.50e-2,                  # P1l
          1.00e-2,                  # P1s
          1.0,                      # P2c
          1.26e-2,                  # P2n
          7.86e-4,                  # P2p
          1.50e-2,                  # P2l
          1.0,                      # P3c
          1.26e-2,                  # P3n
          7.86e-4,                  # P3p
          2.00e-2,                  # P3l
          1.0,                      # P4c
          1.26e-2,                  # P4n
          7.86e-4,                  # P4p
          2.00e-2,                  # P4l
          1.0,                      # Z3c
          1.5e-2,                   # Z3n
          1.67e-3,                  # Z3p
          1.0,                      # Z4c
          1.5e-2,                   # Z4n
          1.67e-3,                  # Z4p
          1.0,                      # Z5c
          1.67e-2,                  # Z5n
          1.85e-3,                  # Z5p
          1.0,                      # Z6c
          1.67e-2,                  # Z6n
          1.85e-3,                  # Z6p
          1.0,                      # R1c
          1.26e-2,                  # R1n
          7.862e-4,                 # R1p
          0.1,                      # R2c
          1.0,                      # R3c
          1.0,                      # R6c
          1.26e-2,                  # R6n
          7.862e-4,                 # R6p
          1.45e-2,                  # R6s
          27060.0,                  # O3c
          2660.0                    # O3h
          ]
    
    # Time span for integration
    t_span = [0, 86400*365*10]
    
    # Set up array to store timer data for each case
    timer_data = {}
    data = []
    
    # Set up aaray to store integration results
    integration_results = []
    
    # Open file and get species removed
    filename = 'reduced_models_species_removed.txt'
    case = 0
    with open(filename) as f:
        for species_removed_str in f:

            # Convert dictionary string to dictionary
            data.append(species_removed_str)
            species_removed = ast.literal_eval(species_removed_str)
            
            # Integrate reduced model and time it
            conc_reduced = copy.copy(c0)
            multiplier = numpy.ones(len(c0))
            for index in species_removed.values():
                conc_reduced[index] = 0.0
                multiplier[index] = 0.0
            start_time = time.time()
            solution_reduced_model = solve_ivp(lambda time, conc: bfm_reduced_rate_eqns(time, conc, multiplier), t_span, conc_reduced, method='RK23')
            end_time = time.time()
            duration = end_time - start_time
            
            # Store simulation duration for each case
            timer_data[case] = duration
            
            # Store integration results
            integration_results.append(solution_reduced_model)
            
            # Increase case number (starts with case 0 for full model)
            case += 1
            
    # Write timer data to output file
    timer_output_file = open("BFM_Reduced_Models_Data/timer_data.txt", "w")
    for case, duration in timer_data.items():
        timer_output_file.write('Case{} simulation time: {} seconds \n'.format(case,duration))
    timer_output_file.close()
    
    # Write integration results to output file
    # Full model
    numpy.savetxt('BFM_Reduced_Models_Data/full_50_t_data.csv', integration_results[0].t)
    numpy.savetxt('BFM_Reduced_Models_Data/full_50_c_data.csv', integration_results[0].y)
    
    # Reduced 44 species model
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_44_t_data.csv', integration_results[1].t)
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_44_c_data.csv', integration_results[1].y)

    # Reduced 21 species model
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_21_t_data.csv', integration_results[2].t)
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_21_c_data.csv', integration_results[2].y)
    
    # Reduced 41 species model
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_41_t_data.csv', integration_results[3].t)
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_41_c_data.csv', integration_results[3].y)
    
    # Reduced 4 species model
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_4_t_data.csv', integration_results[4].t)
    numpy.savetxt('BFM_Reduced_Models_Data/reduced_4_c_data.csv', integration_results[4].y)
