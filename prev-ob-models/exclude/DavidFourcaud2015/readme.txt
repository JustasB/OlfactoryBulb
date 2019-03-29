For any help, contact:
Nicolas Fourcaud-Trocmé: nicolas.fourcaud-trocme@cnrs.fr
François David: francois.david5@free.fr

############################################################

* This model has been run with:

    Python 2.6/2.7 (https://www.python.org/)
    BRIAN 1.4.1 (http://briansimulator.org/)

* Oscillation analysis and time-frequency plots are done with:

    OpenElectrophy 0.2 (http://neuralensemble.org/trac/OpenElectrophy)

    Note that if OpenElectrophy is not installed as a global module,
    OpenElectrophy path can be provided in "oscillation_analysis.py" line 18

###########################################################################

File details:

* model_mitral_clean.py and model_granule_clean.py: 
    contain equations and some fixed parameters of mitral and granule cell models
    Mitral cell model is issued from David et al. Plos Comp Biol (2009),
    Granule cell model is a standard QIF model with f-I curve estimated based on Davison (2001, PhD Thesis)
    
* reseau_mitral_granule_fig_param_dic.py:
    it's the main script to describe full network and launch network simulations. 
    It contains a dictionary of model parameters with all default values and comments on the role of most parameters. 
    The "reseau_mitral_granule" function run the network and return some recordings. It can be safely multiprocessed.
    The "main" section (starting line 229) shows how to run either a single model or a set of models while varying one parameter.
    Each model run can be configured by feeding a distinct dictionay of parameters "param_dict".
    A set of dictionary parameters used in the article is given in "params_for_article_fig.py" and can be easily imported (see commented lines)
    
* oscillation_analysis.py, plot_single_run_from_file.py, plot_multi_run_from_file.py:
    contains functions to analyse gamma/beta oscillations and plot output of network simulations 
    (either detailed output for single simulations, or simplified output as a function of the varied parameters for multiple simulations)
    Note that simulation outputs can be saved in a file (see options in "reseau_mitral_granule_fig_param_dic.py") and 
    later plotted with "plot_single_run_from_file.py" or "plot_multi_run_from_file.py" (see their "main" sections).
    If OpenElectrophy is not installed on the computer, timefrequency plots and oscillation analysis are skipped.
    
* populationstatemonitor.py: is a simple helper function derived from BRIAN simulator "Statemonitor"
