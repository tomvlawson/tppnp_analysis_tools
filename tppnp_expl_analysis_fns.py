"""
tppnp_expl_analysis_fns

This module holds plotting tools to be used with tppnp data, as well as early solar system (ESS) 
plotting tools. Currently only TPPNP, but plans on MPPNP compatibility.
Currently it has only been tested with Sam & Chris' data, but has no issue with data updates as of yet
Plotting can take some time as I have bumped up the resolution of the plots

@author: tomvlawson (Thomas Lawson, Hull)

@TODO:
    Double check with a single file read
        Not too long, but low priotrity
    Ensure compatibility with all tppnp runs
        High priority, give to james??
    create master_table equiv??
        Need to grab data from somewhere
        Create tool to format it automatically?
            IMPORTANT FOR MPPNP DATA SUITE
    Debugging suite:
        Timing
            No point external ising
        printing error
            Debug level??
        key junctions
            Std error, with status update
        rough traj creator shell boundary test
    Rewrite indexing
        No dictionaries
    Label in solarsytem prod factor
    Model ESS comparison plot for many 
    Trajectory Generator from TPPNP data
        Need time, density and temperature profiles
            Not sure if information there with S+C model
        Need initial abundance at a timestep
    Include data aquisition for MPPNP data sets
        Master table generator
        Dictionary creator
        Ensure compatibility with tppnp and mppnp data dictionaries for all the commands
    Temp peak plotter and prog plotter

"""

import sys
import os
# This path will have to be chenged to wherever your tppnp.py is located
sys.path.insert(0, '/Users/thomaslawson/NuGrid/aprNuPPN/frames/tppnp/tools')
import tppnp as t
from nugridpy import utils as u
from nugridpy import nugridse as mp
from nugridpy import ascii_table as at
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import time
from progenitor import progenitor as prog
coll = ['k','b','g','r','y','orange','0.5']*10
symb = ['.','x','+','s','p','h','*','v','>','<']*10

def format_e(n):
    a = '%0.3E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def decay_tau(n_0,del_t,tau):

    return(n_0*np.exp(-(del_t/tau)))

def decay_hl(n_0,del_t,hl):

    return(n_0*np.exp(-((del_t*np.log(2)/hl))))

def ESS_plotter(ESS_path ,iso_in,ref_in,xmin,xmax):
    '''
    ESS_plotter(ESS_path ,iso_in,ref_in,xmin,xmax):no output
    
    This function is used to add ESS ratio data to abundance ratio plots
    NOTE: This function has no inbuilt plt.show() as it is designed to 
        be ran within other functions. If using this on its own add
        the needed plt.show()

    ESS_path:    Path to ESS data, string format
    iso_in:      SLR isotopes, in list for five character format
    ref_in:      reference isotopes, in list for five character format
    xmin:        xmin for plt.hlines, single float needed
    xmax:        xmax for plt.hlines, single float needed
    '''
    ESS_file = open(ESS_path,'r')
    lines = ESS_file.readlines()[2:]
    slr,ref,t12,tau,ess,esse=[],[],[],[],[],[]
    for i in lines:
        slr.append(i.split(',')[0])
        ref.append(i.split(',')[1])
        t12.append(float(i.split(',')[2]))
        tau.append(float(i.split(',')[3]))
        ess.append(float(i.split(',')[4]))
        esse.append(float(i.split(',')[5]))
    ESS_file.close()
    col_count = 0
    for i in range(len(iso_in)):
        iso = iso_in[i]
        ref_iso = ref_in[i]
        for j in range(len(slr)):
            if iso==slr[j]:
                if ref_iso==ref[j]:
                    plt.hlines(ess[j],float(xmin),float(xmax),coll[col_count])
                    plt.hlines(ess[j]+esse[j],float(xmin),float(xmax),coll[col_count],'--')
                    plt.hlines(ess[j]-esse[j],float(xmin),float(xmax),coll[col_count],'--')
                    # print(iso,ref[j],ess[j],esse[j],coll[col_count],j,col_count)
        col_count+=1

def radio_decayed():
    '''
    radio_decayed():radio_decay
    
    This function outputs the radioactive decay paths, indexed by p.radioactives function

    Single output, can be called inline or by asigning variable
    '''
    radio_decay=[\
    ['SI 26', 'P  26', 'S  26', 'CL 26'],\
    [],\
    ['SC 41', 'TI 41', 'V  41', 'CR 41', 'MN 41'],\
    ['V  44', 'CR 44', 'MN 44', 'FE 44', 'CO 44'],\
    ['FE 53', 'CO 53', 'NI 53', 'CU 53', 'ZN 53', 'GA 53'],\
    ['CU 56', 'ZN 56', 'GA 56', 'GE 56'],\
    ['CO 60', 'CU 60', 'ZN 60', 'GA 60', 'GE 60', 'AS 60', 'SE 60'],\
    [],\
    ['RU 97', 'RH 97', 'PD 97', 'AG 97', 'CD 97', 'IN 97', 'SN 97', 'SB 97'],\
    [],\
    ['GA107', 'GE107', 'AS107', 'SE107', 'BR107', 'KR107', 'RB107', 'SR107', 'Y 107', 'ZR107', 'NB107', 'MO107', 'TC107', 'RU107', 'RH107'],\
    ['IN126', 'CD126', 'AG126','PD126','RH126','RU126','TC126','MO126','NB126','ZR126','Y 126','SR126','RB126'],\
    ['TE129', 'SB129', 'SN129', 'IN129', 'CD129','AG129','PD129','RH129','RU129','TC129','MO129','NB129','ZR129','Y 129','SR129'],\
    ['CS135', 'XE135', 'I 135', 'TE135', 'SB135', 'SN135', 'IN135', 'CD135','AG135','PD135','RH135','RU135','TC135','MO135','NB135','ZR135'],\
    ['PM146','EU146', 'GD146', 'TB146', 'DY146', 'HO146', 'ER146', 'GD150', 'TB150', 'DY150', 'HO150', 'ER150', 'DY154', 'HO154', 'ER154','TM154','YB154'],\
    ['LU182', 'YB182', 'TM182', 'ER182', 'HO182', 'DY182', 'TB182', 'GD182','EU182','SM182','PM182','ND182','PR182','CE182','LA182','BA182','CS182','XE182'],\
    ['PB205', 'BI205', 'PO205','AT205','RN205','FR205','RA205','AC205']]
    return(radio_decay)
      
def stable_decayed():
    '''
    radio_decayed():radio_decay
    
    This function outputs the reference decay paths, indexed by p.radioactives function

    Single output, can be called inline or by asigning variable
    '''
    stable_decay=[\
    ['SI 27', 'NE 27', 'NA 27', 'MG 27', 'P  27', 'S  27', 'CL 27', 'AR 27'],\
    ['NE 35', 'NA 35', 'MG 35', 'AL 35', 'SI 35', 'P  35', 'S  35', 'AR 35', 'K  35', 'CA 35', 'SC 35', 'TI 35'],\
    ['SC 40', 'TI 40', 'V  40', 'CR 40', 'MN 40'],\
    ['SC 48', 'V  48', 'CR 48', 'MN 48', 'FE 48', 'CO 48', 'NI 48', 'CU 48'],\
    ['P  55', 'S  55', 'CL 55', 'AR 55', 'K  55', 'CA 55', 'SC 55', 'TI 55', 'V  55', 'CR 55', 'FE 55', 'CO 55', 'NI 55', 'CU 55', 'ZN 55', 'GA 55', 'GE 55'],\
    ['CU 58', 'ZN 58', 'GA 58', 'GE 58', 'AS 58'],\
    ['P  56', 'S  56', 'CL 56', 'AR 56', 'K  56', 'CA 56', 'SC 56', 'TI 56', 'V  56', 'CR 56', 'MN 56', 'CO 56', 'NI 56', 'CU 56', 'ZN 56', 'GA 56', 'GE 56'],\
    ['TC 92', 'RU 92', 'RH 92', 'PD 92', 'AG 92', 'CD 92', 'IN 92'],\
    ['TC 98', 'RH 98', 'PD 98', 'AG 98', 'CD 98', 'IN 98', 'SN 98', 'SB 98'],\
    ['RH 98', 'PD 98', 'AG 98', 'CD 98', 'IN 98', 'SN 98', 'SB 98'],\
    ['GA108', 'GE108', 'AS108', 'SE108', 'BR108', 'KR108', 'RB108', 'SR108', 'Y 108', 'ZR108', 'NB108', 'MO108', 'TC108', 'RU108', 'RH108'],\
    ['IN124', 'CD124'],\
    ['TE127', 'SB127', 'SN127', 'IN127', 'CD127', 'XE127', 'CS127', 'BA127', 'LA127', 'CE127', 'PR127', 'ND127'],\
    ['XE133', 'I 133', 'TE133', 'SB133', 'SN133', 'BA133', 'LA133', 'CE133', 'PR133', 'ND133', 'PM133', 'SM133'],\
    ['EU144', 'GD144', 'TB144', 'DY144', 'HO144', 'GD148', 'TB148', 'DY148', 'HO148', 'ER148', 'TM148'],\
    ['LU180', 'YB180', 'TM180', 'ER180', 'HO180', 'DY180', 'TB180', 'GD180'],\
    ['TL204', 'BI204', 'PO204']]
    return(stable_decay)

def marco_plot_params():
    '''
    Function that sets the plotting style to marco's std settings

    '''

    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 3
    plt.rcParams['xtick.minor.size'] = 5
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 3
    plt.rcParams['ytick.minor.size'] = 5
    plt.rcParams['ytick.minor.width'] = 1
    plt.tick_params(which='major',length=10, width=2)
    plt.tick_params(which='minor',length=4, width=1)

def fivecharformat(stringinput,include_extra_iso=False): # Used when grabbing raw data from iso_massf
    """
    fivecharformat(stringinput,include_extra_iso=False):iso_val
    
    INPUT ARGUMENTS
        stringinput:       String isotope input in form 'El-#'
        include_extra_iso: Attempts to include metastable and ground state isotopes, boolean option
    
    RETURNS
        iso_val:           Isotope formatted from stringinput to five character formatt eg EL123
    """
    import re
    iso_val=[]
    input_string =  (stringinput.upper()).replace("-","") #Remove '-' and Upper case everything
    iso_A_string = re.search(r'\d+', input_string).group() # Find the numerical characters - Atomic Mass
    iso_EL_string = input_string[:input_string.find(iso_A_string[0])] # Find the alphabetical characters
    if len(iso_EL_string) == 2: iso_val = iso_EL_string+iso_A_string.rjust(3) #Create a 5 character length
    elif len(iso_EL_string) == 1: iso_val = iso_EL_string+iso_A_string.rjust(4) # string for the isotope
    elif iso_EL_string.find('*') != -1: # Find metastable isotopes - If wanted
        if include_extra_iso: iso_val = iso_EL_string+iso_A_string.rjust(2)
    elif iso_EL_string.find('TAG') != -1: # Find ground state of Ta - If wanted
        if include_extra_iso: iso_val = iso_EL_string+iso_A_string.rjust(2)
    else: # Anything left - If wanted
        if include_extra_iso: iso_val = iso_EL_string+iso_A_string.rjust(2)
    return (iso_val)

def h5_format(stringinput):
    '''
    h5_format(stringinput):formatted
    
    This returns a five character format element to h5 style: eg 'C-11' or 'Ag-107'
    
    Arguments:
        stringinput: Five character format, string input

    Return: 
        h5 style element name
    '''
    fir=stringinput[0]
    sec=stringinput[1].lower()
    avl=stringinput[2:]
    return(fir+sec+'-'+avl.strip(' '))

def quick_E(n):
    ''' Quick sci not shortening function '''
    return(float('%E' %n))

def var_tool(rat,model_name,norm_slr,model_dict,del_t,ess_dict):
    '''
    var_tool(rat,model_name,norm_slr,model_dict,del_t,ess_dict)

    Tool to create a list of 
    '''
    norm_el = fivecharformat(norm_slr)
    cnt=0
    for j in model_dict[list(model_dict.keys())[0]][0]:
        if j == norm_el:
            decay_norm = decay_tau(model_dict[model_name][1][cnt],del_t,ess_dict[norm_el][2])
            norm_plt = decay_norm/ess_dict[norm_el][5]
            f = ess_dict[norm_el][5]/decay_norm
        cnt+=1
    el = fivecharformat(rat)
    cnt = 0
    for j in model_dict[list(model_dict.keys())[0]][0]:
        if j == el:
            ref_data=decay_tau(model_dict[model_name][1][cnt],del_t,ess_dict[el][2])
            prod=(ref_data/ess_dict[el][5])*f
        cnt+=1
    return(prod)

def elemental_index_dictionary_creator(ref_dict,prog_star_path,file_route):
    '''
    elemental_index_dictionary_creator(ref_dict,prog_star_path,file_route):
                el_index
                
    Creates elemntal indexing dictionary, used in many pltting tools below.
    
    Arguments:
        ref_dict:       A dcitionary to reference, usually a data dictionary
        prog_star_path: Progenitor star path
        file_route:     file path to directory where data is held
    Return:
        el_index: Dictionary of elements na dindex values
            keys: element names in five charamcter format
            values: Index of key element in le index listing, eg  [8] > 'C-11'
    '''
    #Import data from a single file
    p = t.particle_set()
    p.load_final_abundances(file_route+list(ref_dict.keys())[0]+'/output',memmap=True)
    p.load_progenitor(prog_star_path)
    # parts = p.num_particles
    iso_inpu = p.isosfin(1) #generate list of isotopes
    iso_index = []
    for i in range(len(iso_inpu)):
    #     print(iso_inpu[i])
        iso_index.append(str(iso_inpu[i])[2:7])#Strip needed string
    # print(str(iso_inpu[5])[2:7])
    seen=[]
    j=0
    for i in iso_index:#Find meta stable values
        if i in seen:
            el_start = iso_index[j][:2]
            el_end   = iso_index[j][3:]
            iso_index[j] = el_start+'*'+el_end
        else:
            seen.append(i)
        j+=1
    el_index = {iso_index[k]:k for k in range(len(iso_index))} #Create dictionary
    return(el_index)

def tppnp_data_dictionary_creator(file_route,prog_star_path,desired_traj,x_data=True,radio=True,integ=True):
    """
    tppnp_data_dictionary_creator(file_route,prog_star_path,desired_traj,
                                  x_data=True,radio=True):data_dict,radio_dict,
                                  mass_dict,radio_name_index
                                  
    This function creates a set of dictionaries to be used in the plotting 
    tools provided below. It is configured to work with the data provided by
    Sam & Chris, but will likely work with any multiple folder examination.
    
    INPUT ARGUMENTS
    file_route:     File path to directory where tppnp explosion data is held
    prog_star_path: Path to progenitor star .npz file
    desired_traj:   Prefix for the desired trajectory, eg M15s, M20l ...
    x_data=True     Enables abundance data dictionary to be generated
    radio=True      Enables radioactive data dictionary to be generated
    integ=True      Enables integrated dictionary creation
    
    RETURNS
    data_dict:      Abundance data dictionary
        Keys: File names
        Values: particle number, further indexed by isotope abundance. 
            Eg dict[100][8] would provide the abundance of C-11 at particle 100
    radio_dict:     Radioactive data dictionary
        Keys: File names
        Vaues: Indexed as such: Radioactive names, Radioactive abundances, 
                Reference Names, Reference Abundances, Ratio names, 
                ratio abundances, decayed radioactive, decayed reference
    mass_dict: Mass coordinate dictionary
        Keys: File names
        Values: Particle mass coordinate values
    radio_name_index: Radio active name indexing dictionary
        Keys: Ratios, eg 'AL 27/AL 27'
        Values: Index values
    integ_dict: Intgerated abundacnce dictionary
        Keys: File names
        Values: Abundances by element, eg [8] C-11 abundance
    """
    star = prog(file=prog_star_path)
    file_index, radio_ratio,files = [],[],[]
    file_list = os.listdir(file_route)
    for i in file_list:
        if i[:(len(desired_traj))] == desired_traj:
            files.append(file_route+i+'/output')
            file_index.append(i)
    file_data,radio_data,integ_data,radio_names,mass_index=[],[],[],[],[]
    time_tot=0
    # Run through the files provided, creating data sets for each.
    for i in files:
        file_dum=[]
        t_st=time.time()
        p = t.particle_set()
        p.load_final_abundances(i,memmap=True)
        p.load_progenitor(star)
        mass_range = p.m()
        mass_index.append(mass_range)
        ## Based on data required: can add other outputs here potentially
        if x_data: # Abundance Data
            ex_mass_range=mass_range[-1]-mass_range[0]
            for part_num in range(1,p.num_particles+1):
                file_dum.append(p.xisofin(part_num))#/p.num_particles)#/ex_mass_range)
            file_data.append(file_dum)
        if radio: # Radioactive data
            radio_data.append(p.radioactive_yields())
        if integ: # Total output yields
            integ_data.append(p.yields()/ex_mass_range)
        t_en=time.time()
        time_tot+=(t_en-t_st)
        
    print('Average time: %s, Total runs: %s' %(time_tot/len(files),len(files)))
    print('TIME TAKEN: %s mins' %((time_tot)/60.))
    if x_data:
        data_dict = {file_index[i]:file_data[i] for i in range(len(file_index))}
    # Radio requires some fiddling, due to ratio's being required for ease of use
    if integ:
        integ_dict = {file_index[i]:integ_data[i] for i in range(len(file_index))}
    if radio:
        radio_data_t=[]
        j=0
        for i in radio_data:
            radio_ratio.append(i[1]/i[3])
            radio_dum=[]
            for k in range(len(i[0])):
                radio_dum.append(str(i[0][k])+'/'+str(i[2][k]))
            radio_names.append(radio_dum)
            j+=1
        for i in range(len(radio_data)): #Formatting the radio_data for dictionary use
            app_dum=[]
            for j in radio_data[i]:
                app_dum.append(j)
            app_dum.append(radio_names[i])
            app_dum.append(radio_ratio[i])
            if integ: # Here we create decayed values for the radioactive sources
                ell_index = elemental_index_dictionary_creator(data_dict,prog_star_path,file_route)
                s_d = stable_decayed()
                r_d = radio_decayed()
                sum_abu_rad=[]
                for ref in r_d:
                    sum_abu_rad_dum=[]
                    missing=[]
                    for refi in ref:
                        try:
                            sum_abu_rad_dum.append(integ_data[i][ell_index[refi]])
                        except KeyError: missing.append(refi)
                    if len(missing)>0: # Show the missing isotopes form the dacey paths above
                        print('Missing:')
                        print(missing)
                    sum_abu_rad.append(np.sum(sum_abu_rad_dum))
                app_dum.append(sum_abu_rad)
                sum_abu_ref=[]
                for ref in s_d:
                    sum_abu_ref_dum=[]
                    missing=[]
                    for refi in ref:
                        try:
                            sum_abu_ref_dum.append(integ_data[i][ell_index[refi]])
                        except KeyError: missing.append(refi)
                    if len(missing) > 0:
                        print('Missing:')
                        print(missing)
                    sum_abu_ref.append(np.sum(sum_abu_ref_dum))
                app_dum.append(sum_abu_ref)
            radio_data_t.append(app_dum) # rad, xrad, ref, xref, rat, xrat, int_rad, int_ref
        radio_dict = {file_index[i]:radio_data_t[i] for i in range(len(file_index))}
        #radio_ratio_dict = {file_index[i]:radio_ratio[i] for i in range(len(file_index))}   
        radio_name_index = {radio_dict[file_index[0]][4][i]:i for i in range(len(radio_dict[file_index[0]][4]))}
    
    mass_dict = {file_index[i]:mass_index[i] for i in range(len(file_index))}
    
    #TODO: Currently not working if specific dictionaries set to False
    return(data_dict,radio_dict,mass_dict,radio_name_index,integ_dict)
    
def ratio_parameter_plot(radio_dict,radio_name_index,table_dict,radio_plot,plt_x,decay=False,debug=False,\
    ESS_path='',savefig=False):
    '''
    ratio_parameter_plot(data_dict,radio_dict,radio_name_index,table_dict,\
                             radio_plot,plt_x,debug=False)
    
    Plot a parameter plot of abundacnes, El/El_solar. the data is normalised by the mass
    range examined in the tppnp simulation
    
    Arguments:
        data_dict:          Data dictionary for abundance data, generated from
                                tppnp_data_dictionary_creator
        radio_dict:         Radioactives dictionary for abundance data, 
                                generated from tppnp_data_dictionary_creator
        radio_name_index:   Radioactive indexing dictionary for abundance data,
                                generated from tppnp_data_dictionary_creator
        table_dict:         Dictionary that provies simulation parameters,
                                generated by master_table_read
        radio_plot:         Radioactive isotopes that are to be plotted
        plt_x:              Parameter that is to be plotted
        decay=True          Enable decay in abundance plots
        debug=False:        @TODO debugging tools
        
    
    
    '''
    coll = ['k','b','g','r','y','orange','0.5']*10
    symb = ['.','x','+','s','p','h','*','v','>','<']*10
    plt.figure(figsize=(15,7),dpi=400)
    radio_plot_index=[]
    for rad_rat in radio_plot:#Run thorugh the radioactive isotopes to generate ratios
        radio_plot_fcf=fivecharformat(rad_rat)#Format switch
        for i in radio_dict[list(radio_dict.keys())[0]][4]:
            if radio_plot_fcf == i[:5]:#Locate isotope in list of rad/ref
                radio_plot_index.append(radio_name_index[i])

    plot_data_file,plot_data_pltx,plot_data_ratio=[],[],[]
    ## Generate data based on the master table rpovided in the google docs
    for file in radio_dict.keys():
        plot_dum2_file,plot_dum2_pltx,plot_dum2_ratio=[],[],[]
        for ind in radio_plot_index:
            cnt=0
            plot_dum_file,plot_dum_pltx,plot_dum_ratio=[],[],[]
            for i in table_dict['run']:
                if i == file[5:]:
                    plot_dum_file.append(i)
                    plot_dum_pltx.append(float(table_dict[plt_x][cnt]))
                    if decay: 
                        plot_dum_ratio.append(float(radio_dict[file][5][ind])+float(radio_dict[file][7][ind]))
                    else: plot_dum_ratio.append(float(radio_dict[file][5][ind]))
                cnt+=1 
            plot_dum2_file.append(plot_dum_file)
            plot_dum2_pltx.append(plot_dum_pltx)
            plot_dum2_ratio.append(plot_dum_ratio)
        plot_data_file.append(plot_dum2_file)
        plot_data_pltx.append(plot_dum2_pltx)
        plot_data_ratio.append(plot_dum2_ratio)

    plot_data_ratio_sorted,plot_data_pltx_sorted,plot_data_file_sorted=[],[],[]
    for i in range(len(plot_data_ratio[0])):# Reformat data to enable easier plotting | Might be included above
        dum,dumm,dummm=[],[],[]
        for j in range(len(plot_data_ratio)):
            dum.append(plot_data_ratio[j][i][0])
            dumm.append(plot_data_pltx[j][i][0])
            dummm.append(plot_data_file[j][i][0])
        plot_data_ratio_sorted.append(dum)
        plot_data_pltx_sorted.append(dumm)
        plot_data_file_sorted.append(dummm)

    if len(ESS_path)>0:
        ref__ = []
        iso__ = []
        for i in range(len(plot_data_ratio_sorted)):
            iso__.append(radio_dict[list(radio_dict.keys())[0]][4][radio_plot_index[i]][:5])
            ref__.append(radio_dict[list(radio_dict.keys())[0]][4][radio_plot_index[i]][6:])
        # print(iso__)
        # print(ref__)
        # print(ESS_path,radio_plot,'!!!!!!!!!!!',np.min(plot_data_pltx_sorted[0]),np.max(plot_data_pltx_sorted[0]))
        ESS_plotter(ESS_path,iso__,ref__,np.min(plot_data_pltx_sorted[0]),np.max(plot_data_pltx_sorted[0]))

    #Include numpy sorting script

    for i in range(len(plot_data_ratio_sorted)): # Plotting the data sorted above, by ratio defined
        x = [plot_data_pltx_sorted[i][k] for k in np.argsort(plot_data_pltx_sorted[i])]
        y = [plot_data_ratio_sorted[i][k] for k in np.argsort(plot_data_pltx_sorted[i])]
        plt.plot(x,y,color=coll[i],marker=symb[i],\
                 label=radio_dict[list(radio_dict.keys())[0]][4][radio_plot_index[i]])
    # Plotting parameters, std Marco set up
    plt.rcParams['xtick.major.size']  = 10
    plt.rcParams['xtick.major.width'] = 3
    plt.rcParams['xtick.minor.size']  = 5
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size']  = 10
    plt.rcParams['ytick.major.width'] = 3
    plt.rcParams['ytick.minor.size']  = 5
    plt.rcParams['ytick.minor.width'] = 1
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('E$_{exp}$ (10$^{51}$ erg)', fontsize=15)
    plt.ylabel('$Abundance$ $Ratio$', fontsize=15)
    # plt.ylim(10.**(-12.),1.0e0)
    # plt.xlim(min(mass_[0])-.01,max(mass_[0])+.01)
    plt.legend(loc='center left' ,ncol=int(np.ceil(len(plot_data_file_sorted)/15)),\
               bbox_to_anchor=(0.95, 0.5),markerscale=0.8,fontsize=14)
    plt.yscale('log')
    plt.xscale('log')
    plt.tick_params(which='major',length=10, width=2)
    plt.tick_params(which='minor',length=4, width=1)
    plt.gcf().subplots_adjust(bottom=0.1)
    plt.title(str(list(radio_dict.keys())[0][1:3])+' $Solar$ $mass$ $model$')
    name = ''
    for i in radio_plot:
        name+=i
    if savefig:
        plt.savefig(str(list(radio_dict.keys())[0][1:3])+"M_prod_plot"+name+'.png')
    else:
        plt.show()

def radio_file_writer(radio_dict,table_dict,form_str='%7.3E',sep_string=' & ',decay=False):
    '''
    radio_file_writer(data_dict,radio_dict,table_dict,form_str='%7.3E',\
                              sep_string=' & '):
        
    Write a txt file with data from radioactive dictioanries
    
    Arguments:
        radio_dict:    Radioactive dictionary 
        table_dict:    Table dictionary  grenerated from master_table_read
        form_str:      Form that the values will be generated in 
        sep_string:    Seperator symbol between rows
    
    '''
    from nugridpy import ascii_table as att
    info_columns = ['specie','yield [Msun]','specie','yield [Msun]','ratio',r'\\ ']
    i=0
    for file in radio_dict.keys():
        if decay: file_name_table = file+'_radioactives_decay.txt'
        else: file_name_table = file+'_radioactives.txt'
        # headers and format
        jj=0
        if decay:
            for ind in table_dict['run']:
                if ind == file[5:]:
                    headers=[(file+' | Prog Mass: '+table_dict['mprog'][jj]+' | Eexp: '+table_dict['Eexp'][jj]+\
                        ' | mrem: '+table_dict['mrem'][jj]+' | DECAYED')]
                jj+=1
        else:
            for ind in table_dict['run']:
                if ind == file[5:]:
                    headers=[(file+' | Prog Mass: '+table_dict['mprog'][jj]+' | Eexp: '+table_dict['Eexp'][jj]+\
                        ' | mrem: '+table_dict['mrem'][jj])]
                jj+=1
        all_data=[]
        all_data.append(radio_dict[file][0])
        if decay: 
            dec_rad = [(radio_dict[file][6][ii]+radio_dict[file][1][ii]) for ii in range(len(radio_dict[file][6]))]
            all_data.append([form_str % ii for ii in dec_rad])
        else: 
            all_data.append([form_str % ii for ii in radio_dict[file][1]])
        all_data.append(radio_dict[file][2])
        if decay: 
            dec_ref = [(radio_dict[file][7][ii]+radio_dict[file][3][ii]) for ii in range(len(radio_dict[file][6]))]
            all_data.append([form_str % ii for ii in dec_ref])
        else: 
            all_data.append([form_str % ii for ii in radio_dict[file][3]])
        if decay:
            dec_ratio = [(radio_dict[file][6][ii]/radio_dict[file][7][ii]) for ii in range(len(radio_dict[file][6]))]
            all_data.append([form_str % ii for ii in dec_ratio])
        else: 
            all_data.append([form_str % ii for ii in radio_dict[file][5]])
        ### attempt to add the trailing '\\' to a latex table
        final_col=len(radio_dict[file][0])*[r'\\ ']
        all_data.append(final_col)
        att.write(file_name_table,headers,info_columns,all_data,sep=sep_string)
        i+=1
    if decay: print('Files written for %s' %file[:3]+' with decay')
    else: print('Files written for %s' %file[:3])

def master_table_read(path):
    '''
    master_table_read(path):table_dict
    
    Creates a dictionary of the data found within the master table, \
        provided by Sam & Chris in tables
    
    Arguments:
        path:   Path to master table file
    Return:
        table_dict: Dictionary of data from master table
            Keys:   Master table column names
            Values: Data from said columns
    '''
    master_table = open(path,'r')
    lines = master_table.readlines()
    header = lines[0].split()
    lines = lines[1:]
    table_dict = {i:[] for i in header}
    cnt=0
    for i in header:
        for j in lines:
            table_dict[i].append(j.split()[cnt])
        cnt+=1
    return(table_dict)
    
def solar_production_factor_plot(radio_dict,radio_name_index,table_dict,radio_plot,plt_x,solar_path,\
    solar_factor=1.,debug=False,savefig=False,decay=True,prog_path='',ESS_path=''):
    '''
    solar_production_factor_plot(data_dict,radio_dict,radio_name_index,table_dict,
            radio_plot,plt_x,solar_path,solar_factor=1.,debug=False,
            savefig=False,decay=True):
    
    Production factor plotting tool, creates a pordution plot of reference 
    abundance over solar abundance. Currently ONLY for reference production
    factors.
    
    Arguments:
        radio_dict:         Radioactives Dictionary
        radio_name_index:   Radioactives indexing dictionary
        table_dict:         Master table parameters dictionary
        radio_plot:         Radioactive that want to be plotted, list format
        plt_x:              Parameter to be plotted, as found in master table
        solar_path:         Path to solar data
        solar_factor=1.:    Factor for solar data
        debug=False:        Debug mode
        savefig=False:      Save plot produced
        decay=True:         Decay included in abundance
        prog_path:          Progenitor path, by default not turned on
    
    '''
    plt.close() 
    if debug: 
        print('DEBUG MODE')
        t_s = time.time()
    u.solar(solar_path,solar_factor)
    coll = ['k','b','g','r','y','orange','0.5']*10
    symb = ['.','x','+','s','p','h','*','v','>','<']*10
    plt.figure(figsize=(10,5),dpi=400)
    radio_plot_index=[]
    for rad_rat in radio_plot:#Run thorugh the radioactive isotopes to generate ratios
        radio_plot_fcf=fivecharformat(rad_rat)#Format switch
        for i in radio_dict[list(radio_dict.keys())[0]][4]:
            # print(i[6:])
            if radio_plot_fcf == i[6:]:#Locate isotope in list of rad/ref
                radio_plot_index.append(radio_name_index[i])
    # print(radio_plot_index)
    plot_data_file,plot_data_pltx,plot_data_ratio=[],[],[]
    ## Generate data based on the master table rpovided in the google docs
    for file in radio_dict.keys():
        plot_dum2_file,plot_dum2_pltx,plot_dum2_ratio=[],[],[]
        for ind in radio_plot_index:
            cnt=0
            plot_dum_file,plot_dum_pltx,plot_dum_ratio=[],[],[]
            for i in table_dict['run']:
                if i == file[5:]:
                    if table_dict['mbounce'][cnt] != 'NA':
                        plot_dum_file.append(i)
                        plot_dum_pltx.append(float(table_dict[plt_x][cnt]))
                        if decay:
                            plot_dum_ratio.append((float(radio_dict[file][3][ind])+float(radio_dict[file][7][ind]))/\
                                                   (u.solar_abundance[radio_dict[file][2][ind].lower()]))
                                #(float(table_dict['mprog'][cnt])-float(table_dict['mrem'][cnt])))/\
                                
                        else:
                            plot_dum_ratio.append((float(radio_dict[file][3][ind]))/\
                                                   (u.solar_abundance[radio_dict[file][2][ind].lower()]))
                                #(float(table_dict['mprog'][cnt])-float(table_dict['mrem'][cnt])))/\
                                
                cnt+=1 
            plot_dum2_file.append(plot_dum_file)
            plot_dum2_pltx.append(plot_dum_pltx)
            plot_dum2_ratio.append(plot_dum_ratio)
        plot_data_file.append(plot_dum2_file)
        plot_data_pltx.append(plot_dum2_pltx)
        plot_data_ratio.append(plot_dum2_ratio)
    # Reformat data to enable easier plotting | Might be included above
    plot_data_ratio_sorted,plot_data_pltx_sorted,plot_data_file_sorted=[],[],[]
    for i in range(len(plot_data_ratio[0])):
        dum,dumm,dummm=[],[],[]
        for j in range(len(plot_data_ratio)):
            if len(plot_data_ratio[j][i]) != 0:
                dum.append(plot_data_ratio[j][i][0])
                dumm.append(plot_data_pltx[j][i][0])
                dummm.append(plot_data_file[j][i][0])
        plot_data_ratio_sorted.append(dum)
        plot_data_pltx_sorted.append(dumm)
        plot_data_file_sorted.append(dummm)
    #Include numpy sorting script
    # Plotting the data sorted above, by ratio defined
    for i in range(len(plot_data_ratio_sorted)): 
        x = [plot_data_pltx_sorted[i][k] for k in np.argsort(plot_data_pltx_sorted[i])]
        y = [plot_data_ratio_sorted[i][k] for k in np.argsort(plot_data_pltx_sorted[i])]
        plt.plot(x,y,color=coll[i],marker=symb[i],\
                 label=radio_dict[list(radio_dict.keys())[0]][4][radio_plot_index[i]][6:])
    # PROGENITOR
    if len(prog_path)>0:
        prog = np.load(prog_path)
        el_str,seen,prog_prod=[],[],[]
        j=0
        z_progen = prog['z']
        a_progen = prog['a']
        for i in z_progen:
            if j>0: el_str.append(u.get_el_from_z(str(i))+str(a_progen[j]))
            else: el_str.append('NEUT')
            j+=1
        el_list = [fivecharformat(i,include_extra_iso=True) for i in el_str[1:]]
        el_list = [el_str[0]]+el_list
        j=0
        for i in el_list:
            if i in seen:
                el_start = el_list[j][:2]
                el_end   = el_list[j][3:]
                el_list[j] = el_start+'*'+el_end
            else: seen.append(i)
            j+=1
        el_npz_index = {el_list[k]:k for k in range(len(el_list))}
        j=0
        mass_progen = prog['m']
        yps_progen  = prog['yps']
        for i in radio_plot:
            plot_el= fivecharformat(i)
            e_pos = el_npz_index[plot_el]
            x_el = []
            for k in range(len(mass_progen)):
                x_el.append(yps_progen[k][e_pos])
            sol_abu=u.solar_abundance[plot_el.lower()]
            if debug: print(i,'\tInit:'+str(quick_E(np.sum(x_el)/len(mass_progen))),\
                '\tSolar:'+str(quick_E(sol_abu)),'\tProd factor:'+\
                str(quick_E((np.sum(x_el)/len(mass_progen))/sol_abu)))
            plt.plot(np.min(plot_data_pltx_sorted[j]),((np.sum(x_el)/len(mass_progen))/sol_abu),\
                color=coll[j],marker=symb[j],markersize=15.)

            j+=1
    if len(ESS_path)>0:
        ESS_plotter(ESS_path,radio_plot,np.min(plot_data_pltx_sorted[0]),np.min(plot_data_pltx_sorted[0]))
    

    # Eexp_=[0.34,0.54,0.82,1.34]
    # for i in Eexp_:
    #     plt.vlines(i,0,1e5)
    marco_plot_params()
    # Apply correct Units, based on column header
    if plt_x[0] == 'E': units = '$10^{51}$ $ergs$'
    elif plt_x[0] == 't': units = '$s$'
    elif plt_x[0] == 'm': units = '$M\odot$'
    else: units = ' '
    plt.xlabel('$%s$ (%s)' %(plt_x,units), fontsize=15)
    plt.ylabel('$Production$ $Factor$ $El/El_{\odot}$', fontsize=15)
    # plt.ylim(10.**(-12.),1.0e0)
    # plt.xlim(min(mass_[0])-.01,max(mass_[0])+.01)
    plt.legend(loc='center left' ,ncol=int(np.ceil(len(plot_data_file_sorted)/15)),\
               bbox_to_anchor=(1.01, 0.5),markerscale=0.8,fontsize=14)
    plt.yscale('log')
    plt.xscale('log')
    # plt.tick_params(which='major',length=10, width=2)
    # plt.tick_params(which='minor',length=4, width=1)
    plt.gcf().subplots_adjust(bottom=0.1)
    plt.title(str(list(radio_dict.keys())[0][1:3])+' $Solar$ $mass$ $model$')
    # Optional Figure saving
    if savefig:
        plt.savefig(str(list(radio_dict.keys())[0][:3])+'prodfactor.png')
    if debug:
        print(time.time()-t_s)
    plt.show()
    
def mass_co_abuplot(file_read,iso_plot,data_dict,mass_dict,el_index,ESS='',savefig=False,prog_plot='',xlims=[]):
    '''
    mass_co_abuplot(file_read,iso_plot,data_dict,mass_dict,ell_index)

    Plotting tool for mass coordinate abundance plots
    Not currently configured for decay inclusion

    Arguments:
        file_read:  File to plot data from, expectiing single input
        iso_plot:   Isotopes to plot, expecting list format
        data_dict:  Data dictionary, as created by tppnp_data_dictionary_creator
        mass_dict:  mass dictionary, as created by tppnp_data_dictionary_creator
        el_index:   Elemental index, as created by elemental_index_dictionary_creator
        ESS_comparison: Change to a pah to an ESS data file to enable comparison


    '''
    coll = ['k','b','g','r','y','orange','0.5']*10
    plt.figure(figsize=(15,7),dpi=400)
    mass_=[]
    iso__=[]
    coll_cnt = 0
    for iso in iso_plot: #Ordered by particle
        cnt = 0
        mass_dum=[]
        iso__dum=[]
        for i in data_dict[file_read]:
            if iso != 'PROT ':
                iso_fcf = fivecharformat(iso,include_extra_iso=True)
            else:
                iso_fcf = iso
            mass_dum.append(mass_dict[file_read][cnt])
            iso__dum.append(i[el_index[iso_fcf]])
            cnt+=1
        coll_cnt+=1
        mass_.append(mass_dum)
        iso__.append(iso__dum)

    cnt=0
    for i in mass_:
        plt.scatter(i,iso__[cnt],color=coll[cnt],label=iso_plot[cnt],marker='.')
        # print(min(i))
        cnt+=1
    if len(prog_plot)>0:
        m15 = np.load(prog_plot)
        el_str=[]
        j=0
        seen=[]
        for i in m15['z']:
            if j>0:
                el_str.append(u.get_el_from_z(str(i))+str(m15['a'][j]))
            else:
                el_str.append('NEUT')
            j+=1
        el_list = [fivecharformat(i,include_extra_iso=True) for i in el_str[1:]]
        el_list = [el_str[0]]+el_list
        j=0
        for i in el_list:
            if i in seen:
                el_start = el_list[j][:2]
                el_end   = el_list[j][3:]
                el_list[j] = el_start+'*'+el_end
            else:
                seen.append(i)
            j+=1
        el_index = {el_list[k]:k for k in range(len(el_list))}
        j=0
        mass_plot = m15['m']
        yps = m15['yps']
        for i in iso_plot:
            if iso != 'PROT ':
                plot_el = fivecharformat(iso,include_extra_iso=True)
            else:
                plot_el = iso
            # plot_el =  fivecharformat(i,include_extra_iso=False)
            el_pos = el_index[plot_el]
            plot_data=[]
            for k in range(len(mass_plot)):
                plot_data.append(yps[k][el_pos])
            plt.plot(mass_plot,plot_data,color=coll[j],label=i+' preCCSNe',ls='--')
            j+=1

    plt.vlines(min(mass_[0]),10.**(-20.),1.0e0,alpha=0.5,linewidth=5)
    # Plotting parameters, std Marco set up
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 3
    plt.rcParams['xtick.minor.size'] = 5
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 3
    plt.rcParams['ytick.minor.size'] = 5
    plt.rcParams['ytick.minor.width'] = 1
    plt.tick_params(which='major',length=10, width=2)
    plt.tick_params(which='minor',length=4, width=1)
    plt.title(file_read)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('$Initial$ $Mass$ $coordinate$', fontsize=15)
    plt.ylabel('$Mass$ $fraction$', fontsize=15)
    plt.ylim(10.**(-20.),1.0e0)
    # plt.xlim(min(mass_[0])-.01,max(mass_[0])+.01)
    if len(xlims)==2:
        plt.xlim(xlims[0],xlims[1])
        plt.vlines((xlims[0]+xlims[1])/2,0,1)
    else:
        if min(mass_[0])-0.5 > 0:
            # plt.xlim(min(mass_[0])-0.5,max(mass_[0])+.01)
            plt.xlim(1.5,max(mass_[0])+.01)
        else:
            plt.xlim(0,max(mass_[0])+.01)
    plt.legend(loc='center left' ,ncol=int(np.ceil(len(iso_plot)/7)), bbox_to_anchor=(0.95, 0.5),markerscale=0.8,fontsize=14)
    plt.yscale('log')
    plt.gcf().subplots_adjust(bottom=0.1)
    if savefig:
        plt.savefig('abuplot_'+file_read+'.png')
        plt.show()
    else:
        plt.show()
    
def mass_co_ratio_plot(file_read,rad_iso_to_plot,radio_dict,radio_name_index,el_index,data_dict,mass_dict): 
    # Mass plot with radioactive ratios
    ########  INPUTS  ########
    radrat_plot = rad_iso_to_plot
    ##########################
    # pos_index = p.posfin_all()
    # print(np.shape(pos_index))
    
    plt.figure(figsize=(10,5),dpi=400)
    radio_plot_index=[]
    # radio_plot_fcf=fivecharformat(radio_plot)
    for files in file_read:
        radio_plot_dum=[]
        for rad_rat in radrat_plot:#Run thorugh the radioactive isotopes to generate ratios
            radio_plot_fcf=fivecharformat(rad_rat)#Format switch
            for i in radio_dict[files][4]:
                if radio_plot_fcf == i[:5]:#Locate isotope in list of rad/ref
                    radio_plot_dum.append(radio_name_index[i])
        radio_plot_index.append(radio_plot_dum)  #Index position of the isotopes in the list of rad/ref         
    radratio=[]
    cnt=0
    for j in file_read:
        radratdum2=[]
        for i in radio_plot_index[cnt]:
            plot_ratio_dum = radio_dict[j][4][i]
    #         print(el_index[plot_ratio_dum[:5]],'|',el_index[plot_ratio_dum[6:]])
            radratdum=[]
            for k in data_dict[j]:
                radratdum.append(k[el_index[plot_ratio_dum[:5]]]/k[el_index[plot_ratio_dum[6:]]])
            radratdum2.append(radratdum)
        radratio.append(radratdum2)
        cnt+=1
    mass_=[]
    iso__=[]
    file_cnt = 0
    for files in file_read:
        iso_cnt=0
        mass_dum2,iso__dum2=[],[]
        print(np.shape(radio_plot_index))
        for iso in radio_plot_index: #Ordered by particle
            cnt = 0
            mass_dum=[]
            iso__dum=[]
            for i in mass_dict[files]:
        #         iso_fcf = fivecharformat(iso,include_extra_iso=True)
                mass_dum.append(i)
    #             print(cnt,iso_cnt)
                print(file_cnt)
                print(iso_cnt)
                print(np.shape(radratio))
                # print(radratio[file_cnt][iso_cnt])
                iso__dum.append(radratio[file_cnt][iso_cnt])
                cnt+=1
            iso_cnt+=1
            mass_dum2.append(mass_dum)
            iso__dum2.append(iso__dum)
            # print('iso dum %s' %(len(iso__dum[0])))
        mass_.append(mass_dum2)
        iso__.append(iso__dum2)
        file_cnt+=1
    # print(np.shape(mass_),np.shape(iso__))
    
    
    cnt_file,cnt=0,0
    for i in range(len(mass_)):
        cnt_rat=0
        for j in range(len(mass_[i])):
            # print(len(mass_[i][j]),len(iso__[i][j]))
            plt.scatter(mass_[i][j],iso__[i][j][i],color=coll[cnt],\
                        label=file_read[i]+' '+radrat_plot[j],marker='.')
            cnt_rat+=1
            cnt+=1
        cnt_file+=1
        
    # Plotting parameters, std Marco set up
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 3
    plt.rcParams['xtick.minor.size'] = 5
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 3
    plt.rcParams['ytick.minor.size'] = 5
    plt.rcParams['ytick.minor.width'] = 1
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('$Initial$ $Mass$ $coordinate$', fontsize=15)
    plt.ylabel('$Mass$ $ratio$', fontsize=15)
    plt.ylim(1E-5,1e3)
    # plt.xlim(min(mass_[0])-.01,max(mass_[0])+.01)
    plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5))# ,ncol=int(np.ceil(len(iso_plot)/7)),\
    #     bbox_to_anchor=(1.01, 0.5),markerscale=0.8,fontsize=14)
    plt.yscale('log')

    plt.tick_params(which='major',length=10, width=2)
    plt.tick_params(which='minor',length=4, width=1)
    plt.gcf().subplots_adjust(bottom=0.1)
    plt.show()
    
def a_x_prog_comparison(prog_path,integ_dict,model,el_index,x_lim=1e-15,a_lim=210):
    '''
    WRITE DOCSTRING HERE
    '''
    progenitor = np.load(prog_path)
    plt.figure(figsize=(10,5),dpi=400)
    prog_a,prog_yps,prog_m=progenitor['a'],progenitor['yps'],progenitor['m']
    data,a_=[],[]
    for i in list(el_index.keys())[2:]:
        a_.append(float(i[2:].strip(' ').strip('*')))
    for j in range(len(prog_a)):
        data.append(np.sum(prog_yps[:,j])/len(prog_m))
    plt.semilogy(prog_a,data,ls=' ',marker='x',color='k',label='Progenitor')
    plt.semilogy(a_,integ_dict[model][2:],ls=' ',marker='o',color='r',label=model)
    # plt.vlines(92,0,1)
    marco_plot_params()
    plt.ylim(x_lim,1)
    plt.xlim(0,a_lim)
    plt.xlabel('A'),plt.ylabel('Mass fraction')
    plt.legend(loc='center left',bbox_to_anchor=(1.01, 0.5))
    plt.show()

def prog_production_factor_plot(radio_dict,radio_name_index,table_dict,radio_plot,plt_x,prog_path,debug=False,\
    savefig=False,decay=True):
    '''
    solar_production_factor_plot(data_dict,radio_dict,radio_name_index,table_dict,
            radio_plot,plt_x,solar_path,solar_factor=1.,debug=False,
            savefig=False,decay=True):
    
    Production factor plotting tool, creates a pordution plot of reference 
    abundance over solar abundance. Currently ONLY for reference production
    factors.
    
    Arguments:
        radio_dict:         Radioactives Dictionary
        radio_name_index:   Radioactives indexing dictionary
        table_dict:         Master table parameters dictionary
        radio_plot:         Radioactive that want to be plotted, list format
        plt_x:              Parameter to be plotted, as found in master table
        solar_path:         Path to solar data
        solar_factor=1.:    Factor for solar data
        debug=False:        Debug mode
        savefig=False:      Save plot produced
        decay=True:         Decay included in abundance
        prog_path:          Progenitor path, by default not turned on
    
    '''
    if debug: print('DEBUG MODE')
    if len(prog_path)>0:
        prog = np.load(prog_path)
        el_str,seen,prog_prod=[],[],[]
        j=0
        prog_z,prog_a,prog_yps,prog_m=prog['z'],prog['a'],prog['yps'],prog['m']
        for i in prog_z:
            if j>0: el_str.append(u.get_el_from_z(str(i))+str(prog_a[j]))
            else: el_str.append('NEUT')
            j+=1
        el_list = [fivecharformat(i,include_extra_iso=True) for i in el_str[1:]]
        el_list = [el_str[0]]+el_list
        j=0
        for i in el_list:
            if i in seen:
                el_start = el_list[j][:2]
                el_end   = el_list[j][3:]
                el_list[j] = el_start+'*'+el_end
            else: seen.append(i)
            j+=1
        el_npz_index = {el_list[k]:k for k in range(len(el_list))}
        j=0
        for i in radio_plot:
            plot_el= fivecharformat(i)
            e_pos = el_npz_index[plot_el]
            x_el = []
            for k in range(len(prog_m)):
                x_el.append(float(prog_yps[k][e_pos]/len(prog_m)))
            j+=1

    coll = ['k','b','g','r','y','orange','0.5']*10
    symb = ['.','x','+','s','p','h','*','v','>','<']*10
    plt.figure(figsize=(10,5),dpi=400)
    radio_plot_index=[]
    for rad_rat in radio_plot:#Run thorugh the radioactive isotopes to generate ratios
        radio_plot_fcf=fivecharformat(rad_rat)#Format switch
        for i in radio_dict[list(radio_dict.keys())[0]][4]:
            # print(i[6:])
            if radio_plot_fcf == i[6:]:#Locate isotope in list of rad/ref
                radio_plot_index.append(radio_name_index[i])
    # print(radio_plot_index)
    plot_data_file,plot_data_pltx,plot_data_ratio=[],[],[]
    ## Generate data based on the master table rpovided in the google docs
    for file in radio_dict.keys():
        plot_dum2_file,plot_dum2_pltx,plot_dum2_ratio=[],[],[]
        for ind in radio_plot_index:
            cnt=0
            plot_dum_file,plot_dum_pltx,plot_dum_ratio=[],[],[]
            for i in table_dict['run']:
                if i == file[5:]:
                    if table_dict['mbounce'][cnt] != 'NA':
                        plot_dum_file.append(i)
                        plot_dum_pltx.append(float(table_dict[plt_x][cnt]))
                        if decay:
                            plot_dum_ratio.append((float(radio_dict[file][3][ind])+float(radio_dict[file][7][ind]))/\
                                                   (x_el[ind]))
                                #(float(table_dict['mprog'][cnt])-float(table_dict['mrem'][cnt])))/\
                                
                        else:
                            plot_dum_ratio.append((float(radio_dict[file][3][ind]))/\
                                                   (x_el[ind]))
                                #(float(table_dict['mprog'][cnt])-float(table_dict['mrem'][cnt])))/\
                                
                cnt+=1 
            plot_dum2_file.append(plot_dum_file)
            plot_dum2_pltx.append(plot_dum_pltx)
            plot_dum2_ratio.append(plot_dum_ratio)
        plot_data_file.append(plot_dum2_file)
        plot_data_pltx.append(plot_dum2_pltx)
        plot_data_ratio.append(plot_dum2_ratio)
    # Reformat data to enable easier plotting | Might be included above
    plot_data_ratio_sorted,plot_data_pltx_sorted,plot_data_file_sorted=[],[],[]
    for i in range(len(plot_data_ratio[0])):
        dum,dumm,dummm=[],[],[]
        for j in range(len(plot_data_ratio)):
            if len(plot_data_ratio[j][i]) != 0:
                dum.append(plot_data_ratio[j][i][0])
                dumm.append(plot_data_pltx[j][i][0])
                dummm.append(plot_data_file[j][i][0])
        plot_data_ratio_sorted.append(dum)
        plot_data_pltx_sorted.append(dumm)
        plot_data_file_sorted.append(dummm)
    #Include numpy sorting script
    # Plotting the data sorted above, by ratio defined
    for i in range(len(plot_data_ratio_sorted)): 
        x = [plot_data_pltx_sorted[i][k] for k in np.argsort(plot_data_pltx_sorted[i])]
        y = [plot_data_ratio_sorted[i][k] for k in np.argsort(plot_data_pltx_sorted[i])]
        plt.plot(x,y,color=coll[i],marker=symb[i],\
                 label=radio_dict[list(radio_dict.keys())[0]][4][radio_plot_index[i]][6:])
    
    

    # Eexp_=[0.34,0.54,0.82,1.34]
    # for i in Eexp_:
    #     plt.vlines(i,0,1e5)
    marco_plot_params()
    # Apply correct Units, based on column header
    if plt_x[0] == 'E': units = '$10^{51}$ $ergs$'
    elif plt_x[0] == 't': units = '$s$'
    elif plt_x[0] == 'm': units = '$M\odot$'
    else: units = ' '
    plt.xlabel('$%s$ (%s)' %(plt_x,units), fontsize=15)
    plt.ylabel('$Production$ $Factor$ $El/El_{progenitor}$', fontsize=15)
    # plt.ylim(10.**(-12.),1.0e0)
    # plt.xlim(min(mass_[0])-.01,max(mass_[0])+.01)
    plt.legend(loc='center left' ,ncol=int(np.ceil(len(plot_data_file_sorted)/15)),\
               bbox_to_anchor=(1.01, 0.5),markerscale=0.8,fontsize=14)
    plt.yscale('log')
    plt.xscale('log')
    # plt.tick_params(which='major',length=10, width=2)
    # plt.tick_params(which='minor',length=4, width=1)
    plt.gcf().subplots_adjust(bottom=0.1)
    plt.title(str(list(radio_dict.keys())[0][1:3])+' $Solar$ $mass$ $model$')
    # Optional Figure saving
    if savefig:
        plt.savefig(str(list(radio_dict.keys())[0][:3])+'prodfactor.png')
    plt.show()

def ess_dict_creator(file_path):
    '''
    ess_dict_creator(file_path)

    Dictionary creator for ESS data, generated from Lugaro et al 2018. ESS mass fractions generated
        from Lodders et al 2009. If no file is found at file path, a correct file will be created.

    Arguements:
        file_path:  Path to file with ESS data

    Return:
        ess_dict:   Dictionary with ESS data with following indexing format:
                        Keys:   SLR names in five char format
                        Variable index: reference isotope, half life, tau, ESS ratio, ESS error,
                                            ESS mass fraction for SLR, ESS mass fraction for REF
    '''

    try: 
        ESS_file = open(file_path,'r')
    except FileNotFoundError:
        print("NO FILE FOUND AT "+file_path+", CREATING ESS DATA FILE")
        SLR = ['AL 26','BE 10','MN 53','PD107','HF182','I 129','NB 92','NB 92','SM146','CL 36',\
               'FE 60','CA 41','PB205','SN126','CS135','TC 97','TC 97','TC 98','TC 98']

        REF = ['AL 27','BE  9','MN 55','PD108','HF180','I 127','NB 93','MO 92','SM144','CL 35',\
               'FE 56','CA 40','PB204','SN124','CS133','MO 92','RU 98','RU 96','RU 98']
        # Halflife for SM146 taken from N.E.Marks et al 2014
        mass_mod = (0.7112/2.59e10)
        T12 = [0.717,1.388,3.74,6.5,8.9,15.7,34.7,34.7,68,0.301,2.62,0.0994,17.3,0.23,2.3,4.21,\
               4.21,4.2,4.2]
        ESS = [5.23e-5,6e-4,7e-6,6.6e-5,1.018e-4,1.28e-4,1.57e-5,3.2e-5,8.28e-3,2.44e-5,1.01e-8,\
               4.6e-9,1.8e-3,3e-6,2.8e-6,1e-6,1.1e-5,2e-5,6e-5]
        ESSe= [0.13e-5,3e-4,1e-6,0.4e-5,0.043e-4,0.03e-4,0.09e-5,0.3e-5,0.44e-3,0.65e-5,0.27e-8,\
               1.9e-9,1.2e-3,0,0,0,0,0,0]
        TAU = [i/np.log(2) for i in T12]
        N_pre = [8.46e4,0.612,9220,0.359,0.0547,1.1,0.78,0.37,0.008,3920,7.78e5,58500,0.066,\
                0.209,0.371,0.37,0.033,0.099,0.033]
        ESS_mf_ref = [N_pre[i]*int(REF[i][-2:])*mass_mod for i in range(len(REF))]
        ESS_mf_rad = [ESS_mf_ref[i]*ESS[i] for i in range(len(ESS_mf_ref))]
        ##############################
        file = open(file_path,'w')
        file.write('ESS data generated form Lugaro et al 2018, Early solar system mass frac from Lodders 09, CSV filetype \nSLR,REF,T12(Myr),TAU(Myr),ESS,ESSe,ESS mass fractions\n')
        for i in range(len(SLR)):
            line_ = str(SLR[i])+','+str(REF[i])+','+str(T12[i])+','+str(TAU[i])+','+str(ESS[i])+','+str(ESSe[i])+','+str(ESS_mf_rad[i])+','+str(ESS_mf_ref[i])+' \n'
            file.write(line_)

        file.close()
        ESS_file = open(file_path,'r')
    lines = ESS_file.readlines()[2:]
    slr,ref,t12,tau,ess,esse,ess_mf_rad,ess_mf_ref=[],[],[],[],[],[],[],[]
    for i in lines:
        slr.append(i.split(',')[0])
        ref.append(i.split(',')[1])
        t12.append(float(i.split(',')[2])*1e6)
        tau.append(float(i.split(',')[3])*1e6)
        ess.append(float(i.split(',')[4]))
        esse.append(float(i.split(',')[5]))
        ess_mf_rad.append(float(i.split(',')[6]))
        ess_mf_ref.append(float(i.split(',')[7]))
    ESS_file.close()
    ESS_file_index=[slr,ref,t12,tau,ess,esse,ess_mf_rad,ess_mf_ref]
    ess_data=[]
    for i in range(len(slr)):
        dum=[]
        for j in ESS_file_index:
            dum.append(j[i])
        ess_data.append(dum)
    ess_dict={slr[i]:ess_data[i][1:] for i in range(len(slr))}
    return(ess_dict)

def model_ess_slr_plot(norm_slr,plotting_slr,model_name,model_dict,table_dict,ess_dict,del_t=1e6,dpi_in=200,\
    save_fig=False,info=True,change_del_t=False,del_t_list=[],lim=0):
    '''
    model_ess_slr_plot(norm_slr,plotting_slr,model_name,model_dict,table_dict,ess_dict,del_t=1e6,
                            dpi_in=200,save_fig=False,info=True)
    
    Arguements:
        norm_slr:           SLR to normalise the dilution factor by, the base line. Typically Al-26 for CCSNe
        plotting_slr:       SLR to be compared to the normalised SLR
        model_name:         File name for model to be compared to ESS data
        model_dict:         Radioactive dictionary containing model data, as created in this file
        table_dict:         Table dict, as created above
        ess_dict:           ESS dictionary as create in ess_dict_creator
        del_t=1e6:          Time delay for decay calculations
        dpi_in=200:         DPI for plot
        save_fig=False:     Option for saving plotts automatically
        info=True:          Option for model info within plot
        change_del_t=False: Set to true to add multiple delta t values, requires del_t_list
        del_t_list=[]:      List of delta t values to use
    '''
    plt.figure(figsize=(8,5),dpi=dpi_in)
    norm_el = fivecharformat(norm_slr)
    cnt = 0
    for j in table_dict['run']:
        if model_name[5:] == j:
            try:
                mcut = (float(table_dict['mrem'][cnt]))
                Eexp = (float(table_dict['Eexp'][cnt]))
            except ValueError:
                print('No mass cut or explosion energy data in master table')
                mcut,Eexp = 'NA','NA'
        cnt+=1
    if change_del_t: # Potentially change this t a list to allwo for more thatn 2 extra del_t
        del_cont=0
        f_list=[]
        for k in del_t_list:
            cnt=0
            for j in model_dict[list(model_dict.keys())[0]][0]:
                if j == norm_el:
                    decayed_abu = decay_tau(model_dict[model_name][1][cnt],k,ess_dict[norm_el][2])
                    f_list.append(ess_dict[norm_el][5]/decayed_abu)
                cnt+=1
            del_cont+=1
        print(f_list)
        plt.scatter(norm_slr,1,color='k')
    else:
        cnt=0
        for j in model_dict[list(model_dict.keys())[0]][0]:
            if j == norm_el:
                decay_norm = decay_tau(model_dict[model_name][1][cnt],del_t,ess_dict[norm_el][2])
                norm_plt = decay_norm/ess_dict[norm_el][5]
                f = ess_dict[norm_el][5]/decay_norm
            cnt+=1
        plt.scatter(norm_slr,1,color='k')
    if change_del_t: # Potentially change this t a list to allwo for more thatn 2 extra del_t
        ratio_list=[]
        for k in del_t_list:
            ratio_dum=[]
            for rat in plotting_slr:
                el = fivecharformat(rat)
                cnt = 0
                for j in model_dict[list(model_dict.keys())[0]][0]:
                    if j == el:
                        decayed_slr=decay_tau(model_dict[model_name][1][cnt],k,ess_dict[el][2])
                        ratio_dum.append(decayed_slr/ess_dict[el][5])
                    cnt+=1
            ratio_list.append(ratio_dum)
        print(len(ratio_list),len(f_list))

        cnt=0
        for i in ratio_list:
            rat_cnt=0
            for j in i:
                # print(plotting_slr[rat_cnt])
                plt.scatter(plotting_slr[rat_cnt],(j*f_list[cnt]),marker=symb[cnt+1],color=coll[rat_cnt+1])
                rat_cnt+=1
            cnt+=1
    else:
        coll_cnt=1
        for rat in plotting_slr:
            el = fivecharformat(rat)
            cnt = 0
            for j in model_dict[list(model_dict.keys())[0]][0]:
                if j == el:
                    ref_data=decay_tau(model_dict[model_name][1][cnt],del_t,ess_dict[el][2])
                    rat_ratio_plt=(ref_data/ess_dict[el][5])
                cnt+=1
            plt.scatter(rat,(rat_ratio_plt*f),marker=symb[0],s=150,color=coll[coll_cnt])
            coll_cnt+=1
    marco_plot_params()
    plt.yscale('log')
    plt.hlines(1,-0.25,len(plotting_slr)+0.25,linestyles='--')
    plt.xlim(-0.25,len(plotting_slr)+0.25)
    if info:
        plt.scatter(0,0,marker=' ',label='E$_{exp}$  = '+str(Eexp)+' foe')
        plt.scatter(0,0,marker=' ',label='m$_{cut}$  = '+str(mcut)+' M$_{\odot}$')
    plt.title(model_name)
    if change_del_t:
        cnt=0
        for k in f_list:
            plt.scatter(0,0,marker=' ',label='f'+"$_{"+str(del_t_list[cnt]/1e6)+" Myr}$ = "+format_e(k))
            cnt+=1
        cnt=0
        for k in del_t_list:
            plt.scatter(norm_slr,1,marker=symb[cnt+1],color=coll[0],label='$\Delta$t = '+str(k/1e6)+' Myr')
            cnt+=1
        plt.ylabel('Predicted ratio / ESS ratio')
    else:
        plt.scatter(0,0,marker=' ',label='f = '+format_e(f))
        plt.ylabel('Predicted ratio ($\Delta$t='+str(del_t/1e6)+' Myr) / ESS ratio')
    if lim>0:
        plt.ylim(1/lim,lim)
    plt.legend(loc=2,bbox_to_anchor=(0,1),fontsize='xx-small')
    if save_fig:
        plt.savefig(model_name+'_slr_ratio.png')
        print("Saved as "+model_name+"_slr_ratio.png")
        plt.show()
    else:
        plt.show()

def model_ess_slr_data(norm_slr,plotting_slr,model_dict,table_dict,ess_dict,sort_by='MODEL',del_t=1e6):
    '''
    
    '''
    production = []
    if sort_by=='ISO':
        for model_name in model_dict.keys():
            production_dum=[]
            for rat in plotting_slr:
                production_dum.append(var_tool(rat,model_name,norm_slr,model_dict,del_t,ess_dict))
            production.append(production_dum)
    if sort_by=='MODEL':
        for rat in plotting_slr:
            production_dum=[]
            for model_name in model_dict.keys():
                production_dum.append(var_tool(rat,model_name,norm_slr,model_dict,del_t,ess_dict))
            production.append(production_dum)
    return(production)

def tppnp_traj_extractor(part_num,path,file_name,source):
    traj = t.particle_set()
    traj.load_trajectories(path,file_name,source)
    time_traj = traj.time(part_num)
    rho_traj  = traj.dens(part_num)
    temp_traj = traj.temp(part_num)

    at.writeTraj(filename='trajectory_tppnp_test.input',data=[time_traj,temp_traj,rho_traj],\
                ageunit=1, tunit=1, rhounit=0, idNum=0)

def tppnp_trajectory_writer(traj_info,model,mass_coordinate,prog_path,mass_dict,):
    progenitor = np.load(prog_path)

    prog_yps=(progenitor['yps'])
    prog_a=(progenitor['a'])
    prog_z=(progenitor['z'])
    prog_m=(progenitor['m'])

    isos,yps  = [],[]
    for i in range(len(prog_z)):
        if prog_z[i] == 0:
            isos.append(['NEUTRON',1])
        else:
            isos.append([u.get_el_from_z(str(prog_z[i])),str(prog_a[i])])
    #     yps.append(prog_yps[i])
    # Mass coordinate e
    prog_dl = (list(np.absolute(np.array(prog_m)-mass_coordinate)))
    prog_dl_closest = prog_dl.index(min(prog_dl))

    prog_closest = (prog_m[prog_dl_closest])

    model_dl = (list(np.absolute(np.array(mass_dict[model])-mass_coordinate)))
    model_dl_closest = (model_dl.index(min(model_dl)))
    model_closest = (mass_dict[model][model_dl_closest])

    # File creation
    f = open(str(model)+str(mass_coordinate)+'iniabu.dat','w')
    f.write('C Model number | time (yr)  |  Temperature (GK)  | density (cm^-3) \n')
    f.write('C %s  %s  %s  %s \n' %(prog_dl_closest,traj_info['time'][0][prog_dl_closest],traj_info['temp'][0][prog_dl_closest],traj_info['dens'][0][prog_dl_closest]))
    f.write('C Mass coordinate that is really used \n')
    f.write('C %s \n' %(prog_closest))

    for i in range(len(isos))[:-5]:
        f.write('D %s   %s   %s\n' %(isos[i][0].upper(),isos[i][1],prog_yps[prog_dl_closest][i]))
    f.close()



    temp_traj_t9 = np.array(traj_info['temp'])/1e9
    at.writeTraj(filename=str(model)+str(mass_coordinate)+'trajectory.input',data=[traj_info['time'][model_dl_closest],\
        temp_traj_t9[model_dl_closest],traj_info['dens'][model_dl_closest]],ageunit=1, tunit=1, rhounit=0, idNum=0)

def tppnp_trajectory_writer_particle(traj_info,model,particle,prog_path,mass_dict,):
    progenitor = np.load(prog_path)

    prog_yps=(progenitor['yps'])
    prog_a=(progenitor['a'])
    prog_z=(progenitor['z'])
    prog_m=(progenitor['m'])

    isos,yps  = [],[]
    for i in range(len(prog_z)):
        if prog_z[i] == 0:
            isos.append(['NEUTRON',1])
        else:
            isos.append([u.get_el_from_z(str(prog_z[i])),str(prog_a[i])])
    #     yps.append(prog_yps[i])



    mass_coordinate = mass_dict[model][particle]
    # Mass coordinate e
    prog_dl = (list(np.absolute(np.array(prog_m)-mass_coordinate)))
    prog_dl_closest = prog_dl.index(min(prog_dl))
    prog_closest = (prog_m[prog_dl_closest])



    model_dl = (list(np.absolute(np.array(mass_dict[model])-mass_coordinate)))
    model_dl_closest = (model_dl.index(min(model_dl)))
    model_closest = (mass_dict[model][model_dl_closest])
    # print(traj_info['time'][0])




    # File creation
    f = open(str(model)+'_'+str(mass_coordinate)+'_iniabu.dat','w')
    f.write('C Model number | time (yr)  |  Temperature (GK)  | density (cm^-3) \n')
    f.write('C %s  %s  %s  %s \n' %(prog_dl_closest,traj_info['time'][0][prog_dl_closest],traj_info['temp'][0][prog_dl_closest],traj_info['dens'][0][prog_dl_closest]))
    f.write('C Mass coordinate that is really used \n')
    f.write('C %s \n' %(prog_closest))
    for i in range(len(isos))[:-5]:
        f.write('D %s   %s   %s\n' %(isos[i][0].upper(),isos[i][1],prog_yps[prog_dl_closest][i]))
    f.close()




    temp_traj_t9 = np.array(traj_info['temp'])/1e9
    at.writeTraj(filename=str(model)+'_'+str(mass_coordinate)+'_trajectory.input',data=[traj_info['time'][model_dl_closest],\
        temp_traj_t9[model_dl_closest],traj_info['dens'][model_dl_closest]],ageunit=1, tunit=1, rhounit=0, idNum=0)

def temp_peak_ppn_traj_creator(examine_model,examine_mass_co,traj_dict,mass_dict,return_peak=False):
    '''
    This function will generate a file of ppn trajectories and initial abundances from a tppnp
    temperature peak. This peak is taken from the closest particle to the mass coord input

    Arguements:
        examine_model   = Model name that wants to be examined
        examine_mass_co = Mass coordinate of themodel chosen where peak temp will be generated
        traj_dict       = Dictionary of trajectory data, generated in 
        mass_dict       = Dictionary of mass data

    '''
    t_s = time.time()
    model_dl = (list(np.absolute(np.array(mass_dict[examine_model])-examine_mass_co)))
    model_dl_closest = (model_dl.index(min(model_dl)))
    model_closest = (mass_dict[examine_model][model_dl_closest])
    examine_t_peak = max(traj_dict[examine_model[5:]]['temp'][model_dl_closest])
    t_peak=[]
    for i in traj_dict.keys():
        t_dum=[]
        for j in traj_dict[i]['temp']:
            t_dum.append(max(j))
        t_peak.append(t_dum)
    t_peak_closest_particle = []
    for i in t_peak:
        dum = list(np.absolute(np.array(i)-examine_t_peak))
    #     print(min(dum))
        dum = dum.index(min(dum))
    #     print(dum)
        t_peak_closest_particle.append(dum)

    traj_list = list(traj_dict.keys())
    # 
    file_dum = 'ppn_at_'+str(time.time()).split('.')[0]+'_temp'
    print('Files will be located in '+file_dum)
    #MAke a file to work in:
    os.makedirs(file_dum)
    os.chdir(file_dum)
    try:
        #Check the file is clean
        # !rm *.dat
        # !rm *.input
        # !pwd
        cnt =0
        print(os.getcwd())
        for i in t_peak_closest_particle:
            tppnp_trajectory_writer_particle(traj_dict[traj_list[cnt]],'M15s_'+traj_list[cnt],i,\
                    "../npz/M15s.20180607.npz",mass_dict)   
            cnt+=1
        print('Files in: '+file_dum)
        os.chdir('..')
        print(str(time.time()-t_s)+'s')
    except:
        raise Exception('SOMETHING WENT WRONG, for bugfixing make sure to comment out the os.chdir() command')
        # os.chdir() is found before try statement
        os.chdir('..')
    if return_peak:
        return(examine_t_peak)

def dens_peak_ppn_traj_creator(examine_model,examine_mass_co,traj_dict,mass_dict):
    '''
    This function will generate a file of ppn trajectories and initial abundances from a tppnp
    density peak. This peak is taken from the closest particle to the mass coord input

    Arguements:
        examine_model   = Model name that wants to be examined
        examine_mass_co = Mass coordinate of themodel chosen where peak temp will be generated
        traj_dict       = Dictionary of trajectory data, generated in 
        mass_dict       = Dictionary of mass data

    '''
    t_s = time.time()
    model_dl = (list(np.absolute(np.array(mass_dict[examine_model])-examine_mass_co)))
    model_dl_closest = (model_dl.index(min(model_dl)))
    model_closest = (mass_dict[examine_model][model_dl_closest])
    examine_t_peak = max(traj_dict[examine_model[5:]]['dens'][model_dl_closest])
    t_peak=[]
    for i in traj_dict.keys():
        t_dum=[]
        for j in traj_dict[i]['dens']:
            t_dum.append(max(j))
        t_peak.append(t_dum)
    t_peak_closest_particle = []
    for i in t_peak:
        dum = list(np.absolute(np.array(i)-examine_t_peak))
    #     print(min(dum))
        dum = dum.index(min(dum))
    #     print(dum)
        t_peak_closest_particle.append(dum)

    traj_list = list(traj_dict.keys())
    # 
    file_dum = 'ppn_at_'+str(time.time()).split('.')[0]+'_density'
    print('Files will be located in '+file_dum)
    #MAke a file to work in:
    os.makedirs(file_dum)
    os.chdir(file_dum)
    try:
        #Check the file is clean
        # !rm *.dat
        # !rm *.input
        # !pwd
        cnt =0
        print(os.getcwd())
        for i in t_peak_closest_particle:
            tppnp_trajectory_writer_particle(traj_dict[traj_list[cnt]],'M15s_'+traj_list[cnt],i,\
                    "../npz/M15s.20180607.npz",mass_dict)   
            cnt+=1
        print('Files in: '+file_dum)
        os.chdir('..')
        print(time.time()-t_s)
    except:
        raise Exception('SOMETHING WENT WRONG, for bugfixing make sure to comment out the os.chdir() command')
        # os.chdir() is found before try statement
        os.chdir('..')

def trajectory_dictionary_creator(path_to_traj,desired,model_type):
    file_list = os.listdir(path_to_traj)
    files,file_traj_index,dat_files=[],[],[]
    for i in file_list:
        if i[:(len(desired))] == desired:
            files.append(path_to_traj+i)
            dat_files.append('traj'+desired[-2:]+str(i[(len(desired)):])+'.dat')
            file_traj_index.append(i)
    data_traj,error=[],[]
    t_s = time.time()
    keys = ['time','temp','dens']
    for i in range(len(files)):
        p = t.particle_set()
        p.load_trajectories(files[i],dat_files[i],model_type)
        rho_part_traj,temp_part_traj,time_part_traj=[],[],[]
        for j in range(1,p.num_particles+1):
            rho_part_traj.append(p.dens(j))
            temp_part_traj.append(p.temp(j))
            time_part_traj.append(p.time(j))
    #     data_traj.append([time_part_traj,temp_part_traj,rho_part_traj])
        if len(rho_part_traj) >0:
            data_listed = [time_part_traj,temp_part_traj,rho_part_traj]
            data_traj_dict = {keys[k]:data_listed[k] for k in range(len(data_listed))}
            data_traj.append(data_traj_dict)
            print(str(i+1)+'/'+str(len(files)))
        else:
            data_listed = [[0,0],[0,0],[0,0]]
            data_traj_dict = {keys[k]:data_listed[k] for k in range(len(data_listed))}
            data_traj.append(data_traj_dict)
            error.append(i)
            print(str(i+1)+'/'+str(len(files)))
    # print(np.shape(data))
    print(len(file_traj_index))
    print('---')
    print(len(data_traj))
    traj_dict = {file_traj_index[i]:data_traj[i] for i in range(len(file_traj_index))}
    print('Issues with following files:')
    for i in error:
        print(file_list[i])

    print('Time taken: %s min' %((time.time()-t_s)/60.))
    return(traj_dict)

def peak_ppn_traj_creator(examine_peak,traj_dict,mass_dict,prog_path,peak_type='temp'):
    '''
    This function will generate a file of ppn trajectories and initial abundances from a tppnp
    temperature peak. This peak is taken from the closest particle to the mass coord input

    Arguements:
        examine_model   = Model name that wants to be examined
        examine_mass_co = Mass coordinate of themodel chosen where peak temp will be generated
        traj_dict       = Dictionary of trajectory data, generated in 
        mass_dict       = Dictionary of mass data

    '''
    t_s = time.time()
    # model_dl = (list(np.absolute(np.array(mass_dict[examine_model])-examine_mass_co)))
    # model_dl_closest = (model_dl.index(min(model_dl)))
    # model_closest = (mass_dict[examine_model][model_dl_closest])
    # examine_t_peak = max(traj_dict[examine_model[5:]]['temp'][model_dl_closest])
    peak=[]
    traj_list = list(traj_dict.keys())
    # print(traj_dict.keys())
    for i in traj_dict.keys():
        t_dum=[]
        # print(traj_dict[i]['temp'])
        for j in traj_dict[i]['temp']:
            t_dum.append(max(j))
        peak.append(t_dum)
    peak_closest_particle = []
    # print(peak)
    for i in peak:
        try:
            dum = list(np.absolute(np.array(i)-examine_peak))
            dum = dum.index(min(dum))
        except ValueError:
            dum = 1
            print('MISSING DATA IN FILE')
        # print(dum)
        peak_closest_particle.append(dum)
        


    
    # 
    file_dum = 'ppn_at_'+str(time.time()).split('.')[0]+'_temp'
    print('Files will be located in '+file_dum)
    #MAke a file to work in:
    os.makedirs(file_dum)
    os.chdir(file_dum)
    #Check the file is clean
    # !rm *.dat
    # !rm *.input
    # !pwd
    print(os.getcwd())
    try:
        cnt =0
        for i in peak_closest_particle:
            tppnp_trajectory_writer_particle(traj_dict[traj_list[cnt]],list(mass_dict.keys())[0][:5]+traj_list[cnt]\
                ,i,"../"+prog_path,mass_dict)
            cnt+=1
    except:
        os.chdir('..')
        raise Exception('SOMETHING WENT WRONG, please check current work directory is correct with os.getcwd()')
        # print('missing file')
    print('Files in: '+file_dum)
    os.chdir('..')
    print(time.time()-t_s)

    # os.chdir() is found before try statement
    





