#source_generator.py>

import openmc

def create_source(source_dict):

    """Creates a list of openmc.Source() specified from a yaml file

    Inputs:
    -------
    source_dict: dict, a dictionary containing data dictionaries


    Returns:
    -------
    sources_out: list, a list of openmc.Source()
    """

    sources_out = [] # Initialize sources list
    
    SORT_ORDER = {"a": 0, "b": 1, "x": 2, "y": 3, "r": 4,\
                  "theta": 5, "p": 6, "interpolation": 7,\
                  "coefficients": 8, "probability": 9,\
                  "distribution": 10, "mean_value": 11,\
                  "std_dev": 12, "e0": 13, "m_rat": 14,\
                   "kt": 15, "mu": 16, "phi": 17, "z": 18,\
                  "reference_uvw": 19, "origin": 20,\
                  "lower_left": 21, "upper_right": 22,\
                  "only_fissionable": 23, "xyz": 24,\
                  "origin": 25, "reference_uvw": 26} 
                     
    for source_name, distributions in source_dict.items():

        source = openmc.Source()
       
       # The distributions that define the source will be
       # separated out by type and iterated over, adding
       # them to the source definition.

        dist_name_list = ['space', 'angle', 'energy']
        dist_list = [distributions['space'],\
                     distributions['angle'],\
                     distributions['energy']]
        
        for dist in dist_list: 
             
            dist_variables = [variables for variables in dist.keys()\
                              if variables not in ['distribution']]
            dist_variables = sorted(dist_variables,\
                                    key=lambda val:SORT_ORDER[val]) 
                 
            # Create a structure for each component P.D.F. to
            # pass into the top-level distribution.

            pdf_data = {} 
            dist_args = {}
        
            for variable in dist_variables:
            
                # Each distribution (energy, angle, space) can either
                # be made up of distributions or be a single
                # distribution itself, which must be considered
                # when generating the instance of the distribution
                
                # Check if the variables the distribution are
                # themselves distributions
               
                if isinstance(dist[variable], (dict, list))\
                   and 'distribution' in dist[variable]:

                    pdf_func = getattr(openmc.stats,
                                       dist[variable]['distribution'])
                               
                    del dist[variable]['distribution']
            
                    # Now need to iterate over the data in this variable
            
                    # Each distribution in OpenMC needs to be it's own
                    # instance of a distribution class. Here each of
                    # the distributions is individually instantiated.
            
                    pdf_data = dist[variable]

                    pdf_func = pdf_func(**pdf_data)
                     
                    dist_args[variable] = pdf_func
                   
                else: # otherwise assign the variable as an argument
                   
                   dist_args[variable] = dist[variable]
            
            # Lastly, instantiate the top-level distribution using
            # the distributions and/or variables that define it,
            # generated previously
                    
            dist_func = getattr(openmc.stats, dist['distribution'])
            dist_func = dist_func(**dist_args)
            
            setattr(source, dist_name_list[dist_list.index(dist)], dist_func)
      
            # Assign these classes and the rest of the data to
            # the appropriate variables in the source
        
            if 'strength' in distributions.keys():
            
                strength = float(distributions['strength'])
                
                source.strength = strength
        
            if 'particle' in distributions.keys():
        
                source.particle = distributions['particle']
                
        sources_out.append(source)
        
    return(sources_out)
