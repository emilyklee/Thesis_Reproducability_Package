import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import pkg_resources
import cantera as ct
import logging
import yaml

from pymars import pymars
from pymars import drgep
from pymars import sampling
from pymars import reduce_model
from pymars import soln2cti
    
    
def relative_location(file):
    file_path = os.path.join(file)
    return pkg_resources.resource_filename(__name__, file_path)

def calc_percent_error(new_array, old_array):
    """ Calculates the percent error between two matricies.
    This is used for the calculating the New Method's direct interaction coeffs
    """
    
    percent_error = []
    if len(new_array)==len(old_array):
        for i in range(len(new_array)):
            if new_array[i] == old_array[i]:
                percent_error.append(0.0)
            else:
                percent_error.append(abs(new_array[i] - old_array[i])/abs(old_array[i])*100)
            
    return percent_error

def get_drgep_oic(model_file, species_targets, ignition_conditions, psr_conditions=[], 
                  flame_conditions=[], phase_name='', num_threads=1, path=''):
    """ Gets the overall interaction coeffs for the DRGEP method """
    
    solution_object = ct.Solution(model_file)
    
    # first, sample thermochemical data and generate metrics for measuring error
    # (e.g, ignition delays). Also produce adjacency matrices for graphs, which
    # will be used to produce graphs for any threshold value.
    sampled_metrics, sampled_data = sampling.sample(model_file, ignition_conditions, psr_conditions, 
                                                    flame_conditions, phase_name, num_threads, path)
    
    matrices = []
    for state in sampled_data:
        matrices.append(drgep.create_drgep_matrix((state[0], state[1], state[2:]), solution_object))
    
    # For DRGEP, find the overall interaction coefficients for all species 
    # using the maximum over all the sampled states
    importance_coeffs = drgep.get_importance_coeffs(
        solution_object.species_names, species_targets, matrices
        )
    
    return importance_coeffs


def get_new_method_oic(model_file, species_targets, ignition_conditions, psr_conditions=[], 
                  flame_conditions=[], phase_name='', num_threads=1, path=''):
    """ Gets the overall interaction coeffs for the DRGEP method """
    
    solution_object = ct.Solution(model_file)
    
    # first, sample thermochemical data and generate metrics for measuring error
    # (e.g, ignition delays). Also produce adjacency matrices for graphs, which
    # will be used to produce graphs for any threshold value.
    sampled_metrics, sampled_data = sampling.sample(model_file, ignition_conditions, psr_conditions, 
                                                    flame_conditions, phase_name, num_threads, path)
    
    matrices = []
    for state in sampled_data:
        matrices.append(get_new_method_dic_matrix((state[0], state[1], state[2:]), solution_object))
    
    # For DRGEP, find the overall interaction coefficients for all species 
    # using the maximum over all the sampled states
    importance_coeffs = drgep.get_importance_coeffs(
        solution_object.species_names, species_targets, matrices
        )
    
    return importance_coeffs


def get_new_method_dic_matrix(state, solution):
    """ Calculates the direct interaction coeffs for the New Method """
    
    temp, pressure, mass_fractions = state
    solution.TPY = temp, pressure, mass_fractions
    
    
    # Get number of species and reactions
    num_species = solution.n_total_species
    num_reactions = solution.n_reactions
        
    # Net production rate (creation - destruction)
    net_production_rates_og = solution.net_production_rates
    
    # Stoichiometric coeffs
    reactant_stoich_coeffs = solution.reactant_stoich_coeffs()
    product_stoich_coeffs = solution.product_stoich_coeffs()
    
    # Find new method direct interaction coeffs
    percent_error_matrix = np.zeros([num_species,num_species])
    
    for j in range(num_species):
        
        # Remove species j
        # Identify all rxns containing species j and store rxn index to a list
        rxn_list = []
        for rxn in range(num_reactions):
            if reactant_stoich_coeffs[j,rxn]>0:
                rxn_list.append(rxn)
            elif product_stoich_coeffs[j,rxn]>0:
                rxn_list.append(rxn)
        # print(j,rxn_list)
        
        # For each reaction that contains species j, use set_multiplier() to "disable" the rxn
        for rxn in rxn_list:
            solution.set_multiplier(0,rxn)
            
        # Get new net production rates
        net_production_rates_new = solution.net_production_rates
        
        # Calculate percent error in orginal versus new net production rate for each species
        # and store values in a matrix
        p_error = calc_percent_error(net_production_rates_new, net_production_rates_og)
        percent_error_matrix[:,j]=p_error
        
        # set multiplier back to one
        for rxn in rxn_list:
            solution.set_multiplier(1,rxn)
        
    # Find row max to normalize matrix
    row_max = np.amax(percent_error_matrix, axis=1)
    
    # Normalize percent_error_matrix by column max which is new_dic_matrix
    new_dic_matrix_row = np.zeros([num_species, num_species])
    for i in range(num_species):
        if row_max[i] == 0:
            new_dic_matrix_row[i,:] = 0.0
        else:
            new_dic_matrix_row[i,:] = percent_error_matrix[i,:]/row_max[i]
            
    # Set diagonals to zero to avoid self_directing graph edges
    np.fill_diagonal(new_dic_matrix_row, 0.0)
    
    return new_dic_matrix_row


    
def reduce_drgep(model_file, species_safe, threshold, importance_coeffs, ignition_conditions, 
                 sampled_metrics, phase_name='', previous_model=None, num_threads=1, path=''
                 ):
    """Given a threshold and DRGEP coefficients, reduce the model and determine the error.

    Parameters
    ----------
    model_file : str
        Filename for model being reduced
    species_safe : list of str
        List of species to always be retained
    threshold : float
        DRG threshold for trimming graph
    importance_coeffs : dict
        Dictionary with species and their overall interaction coefficients.
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    sampled_metrics: numpy.ndarray
        Global metrics from original model used to evaluate error
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    previous_model : ReducedModel, optional
        Model produced at previous threshold level; used to avoid repeated work.
    num_threads : int, optional
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    path : str, optional
        Optional path for writing files

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

    """
    solution = ct.Solution(model_file, phase_name)
    species_removed = [sp for sp in solution.species_names
                       if importance_coeffs[sp] < threshold 
                       and sp not in species_safe
                       ]
    
    if (previous_model and 
        len(species_removed) == solution.n_species - previous_model.model.n_species
        ):
        return previous_model

    # Cut the exclusion list from the model.
    reduced_model = reduce_model.trim(
        model_file, species_removed, f'reduced_{model_file}', phase_name=phase_name
        )
    reduced_model_filename = soln2cti.write(
        reduced_model, f'reduced_{reduced_model.n_species}.cti', path=path
        )

    reduced_model_metrics = sampling.sample_metrics(
        reduced_model_filename, ignition_conditions, phase_name=phase_name, 
        num_threads=num_threads, path=path
        )
    error = sampling.calculate_error(sampled_metrics, reduced_model_metrics)

    return reduce_model.ReducedModel(
        model=reduced_model, filename=reduced_model_filename, error=error
        )


def run_drgep(model_file, importance_coeffs, ignition_conditions, psr_conditions, flame_conditions, 
              error_limit, species_targets, species_safe, phase_name='',
              threshold_upper=None, num_threads=1, path=''
              ):
    """Main function for running DRGEP reduction.
    
    Parameters
    ----------
    model_file : str
        Original model file
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    psr_conditions : list of InputPSR
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame
        List of laminar flame simulation conditions.
    error_limit : float
        Maximum allowable error level for reduced model
    species_targets : list of str
        List of target species names
    species_safe : list of str
        List of species names to always be retained
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    threshold_upper : float, optional
        Upper threshold (epsilon^*) to identify limbo species for sensitivity analysis
    num_threads : int, optional
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    path : str, optional
        Optional path for writing files

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata
    
    """
    solution = ct.Solution(model_file, phase_name)
    assert species_targets, 'Need to specify at least one target species.'

    # first, sample thermochemical data and generate metrics for measuring error
    # (e.g, ignition delays). Also produce adjacency matrices for graphs, which
    # will be used to produce graphs for any threshold value.
    sampled_metrics, sampled_data = drgep.sample(
        model_file, ignition_conditions, phase_name=phase_name, num_threads=num_threads, path=path
        )
    
    matrices = []
    for state in sampled_data:
        matrices.append(drgep.create_drgep_matrix((state[0], state[1], state[2:]), solution))

    # For DRGEP, find the overall interaction coefficients for all species 
    # using the maximum over all the sampled states
    # importance_coeffs = drgep.get_importance_coeffs(
    #     solution.species_names, species_targets, matrices
    #     )

    # begin reduction iterations
    logging.info('Beginning reduction loop')
    logging.info(45 * '-')
    logging.info('Threshold | Number of species | Max error (%)')
    
    # make lists to store data
    threshold_data = []
    num_species_data = []
    error_data = []

    # start with detailed (starting) model
    previous_model = reduce_model.ReducedModel(model=solution, filename=model_file, error=0.0)

    first = True
    error_current = 0.0
    threshold = 0.01
    threshold_increment = 0.01
    while error_current <= error_limit:
        reduced_model = reduce_drgep(
            model_file, species_safe, threshold, importance_coeffs, ignition_conditions, 
            sampled_metrics, phase_name=phase_name, previous_model=previous_model, 
            num_threads=num_threads, path=path
            )
        error_current = reduced_model.error
        num_species = reduced_model.model.n_species

        # reduce threshold if past error limit on first iteration
        if first and error_current > error_limit:
            error_current = 0.0
            threshold /= 10
            threshold_increment /= 10
            if threshold <= 1e-6:
                raise SystemExit(
                    'Threshold value dropped below 1e-6 without producing viable reduced model'
                    )
            logging.info('Threshold value too high, reducing by factor of 10')
            continue
        
        logging.info(f'{threshold:^9.2e} | {num_species:^17} | {error_current:^.2f}')
        
        # store data
        threshold_data.append(threshold)
        num_species_data.append(num_species)
        error_data.append(error_current)

        threshold += threshold_increment
        first = False

        # cleanup files
        if previous_model.model.n_species != reduced_model.model.n_species:
            os.remove(reduced_model.filename)
        
        previous_model = reduce_model.ReducedModel(
            model=reduced_model.model, filename=reduced_model.filename, 
            error=reduced_model.error, limbo_species=reduced_model.limbo_species
            )
    
    if reduced_model.error > error_limit:
        threshold -= (2 * threshold_increment)
        reduced_model = reduce_drgep(
            model_file, species_safe, threshold, importance_coeffs, ignition_conditions, 
            sampled_metrics, phase_name=phase_name, num_threads=num_threads, path=path
            )
    else:
        soln2cti.write(reduced_model, f'reduced_{reduced_model.model.n_species}.cti', path=path)

    if threshold_upper:
        for sp in reduced_model.model.species_names:
            if importance_coeffs[sp] < threshold_upper and (sp not in species_safe):
                reduced_model.limbo_species.append(sp)
    
    logging.info(45 * '-')
    logging.info('Reduction complete.')
    logging.info(f'Skeletal model: {reduced_model.model.n_species} species and '
                 f'{reduced_model.model.n_reactions} reactions.'
                 )
    logging.info(f'Maximum error: {reduced_model.error:.2f}%')
    logging.info('Final reduced model saved as ' + reduced_model.filename)
    
    return threshold_data, num_species_data, error_data


if __name__ == '__main__':

    # General information for both methods
    # input .yaml filename
    yaml_file = 'ep46_mech_drgep_nosa.yaml'
    
    # input dict
    with open(yaml_file, 'r') as the_file:
        input_dict = yaml.safe_load(the_file)
    # Parse dictionary and check for consistency and correctness
    reduction_inputs = pymars.parse_inputs(input_dict)
    # Get all inputs
    ignition_conditions = reduction_inputs.ignition_conditions
    model_file = reduction_inputs.model
    error_limit = reduction_inputs.error
    psr_conditions = reduction_inputs.psr_conditions
    flame_conditions = reduction_inputs.flame_conditions
    species_targets = reduction_inputs.target_species
    species_safe = reduction_inputs.safe_species
    upper_threshold = reduction_inputs.upper_threshold
    phase_name = reduction_inputs.phase_name
    
    # other info not from input file
    num_threads = 1
    path = ''
    
    #------------------------------ DRGEP ------------------------------------
    logging.info('DRGEP Method:')
    # Get DRGEP overall interactiopn coefficients
    drgep_overall_interaction_coeffs = get_drgep_oic(
        model_file, species_targets, ignition_conditions, psr_conditions, flame_conditions, phase_name, num_threads, path
        )
    
    threshold_data_drgep, num_species_data_drgep, error_data_drgep = run_drgep(model_file, drgep_overall_interaction_coeffs, ignition_conditions, psr_conditions, flame_conditions, error_limit, species_targets, species_safe)
    
    #---------------------------- New Method ---------------------------------
    logging.info('Modified DRGEP Method:')
    upper_threshold_new_method = 1e-7

    # Get the overall interaction coeffs (capable of multiple ignition conditions)
    new_method_overall_interaction_coeffs = get_new_method_oic(model_file, species_targets, ignition_conditions, psr_conditions, flame_conditions, phase_name, num_threads, path)

    # Run reduction on new method
    threshold_data_new_method, num_species_data_new_method, error_data_new_method = run_drgep(model_file, new_method_overall_interaction_coeffs, ignition_conditions, psr_conditions, flame_conditions, error_limit, species_targets, species_safe, phase_name, upper_threshold_new_method, num_threads, path)
