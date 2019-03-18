matlab code for
Chow S-F, Wick SD, Riecke H (2012) Neurogenesis Drives Stimulus
Decorrelation in a Model of the Olfactory Bulb. PLoS Comput Biol 8(3):
e1002398. doi:10.1371/journal.pcbi.1002398

run main.m for figures similar to fig.2, but with only 26 channels
    set sim = 10 for 442 channels

run main2.m for figures similar to fig.9
	enrich = [1 2]; for related enrichment
    enrich = [3 4]; for unrelated enrichment

figures:

    1: input patterns and correlations
    2: output patterns and correlations, connection
    3: GC information
    4: correlation
        blue: mean correlation of all paris
        red: mean correlation of tracked pairs (see parameter tracking)
    101: 2d representation of input patterns
    102: 2d representation of output patterns

parameters:

    sim: level of down sampling
        for the default input set, sim 40 -> 26 channels, sim 10 ->
        442 channels
    odor_names: file names of patterns from http://gara.bio.uci.edu/
    choose: odors chosen to be in the training set

    non_lin: 0 for linear network, 1 for rectified nonlinear network
             (much slower)
    conn: number of connection each GC makes (mean # connection for
          prob_conn = 1)
    CS: coupling strength

    ts: threshold of survival function
    gamma: steepness of survival function
    th: minimal GC activity that would count towards survival
    rm, rg: thresholds for rectifer, only works for non_lin == 1

    cont_density: 0 for discrete GC population, 1 for population description
        population description version runs very slow for large network
    exp_time: total experiment time
    step: plotting/output interval
    dt: for equation stepping

    tracking:
    track mean correlation for a subset of the odor pairs

    We gratefully acknowledge the support of NSF grant DMS-0719944

20121127 matlab code in main.m main2.m modified by replacing ~ with
the variable "ignore" for backwards compatibility with matlab versions
before R2009b.
