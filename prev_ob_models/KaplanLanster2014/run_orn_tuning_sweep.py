import os
import numpy as np


def run_full_sweep():

    gor_min_min_sweep = 5e-5
    gor_min_max_sweep = 5e-5
    gor_max_min_sweep = .7e-3
    gor_max_max_sweep = 1.1e-3

    gkcag_min_min_sweep = 5e-3
    gkcag_min_max_sweep = 5e-3
    gkcag_max_min_sweep = 5e-2
    gkcag_max_max_sweep = 5e-2

    gcal_min_min_sweep = 3e-5
    gcal_min_max_sweep = 3e-5
    gcal_max_min_sweep = 0.8e-5
    gcal_max_max_sweep = 0.8e-5

    gleak_min_min_sweep = 8e-5
    gleak_min_max_sweep = 8e-5
    gleak_max_min_sweep = 1.2e-4
    gleak_max_max_sweep = 1.2e-4

    #gor_max_sweep = 1e-3
    #gkcag_min_sweep = 1e-3
    #gkcag_max_sweep = 1e-2
    #gcal_min_sweep = 1e-4
    #gcal_max_sweep = 1e-3
    #gleak_min_sweep = 1e-4
    #gleak_max_sweep = 1e-5

    n_gor_steps = 1
    gor_min_range = np.linspace(gor_min_min_sweep, gor_min_max_sweep, n_gor_steps)
    n_gor_steps = 5
    gor_max_range = np.linspace(gor_max_min_sweep, gor_max_max_sweep, n_gor_steps)
    n_gkcag_steps = 1
    gkcag_min_range = np.linspace(gkcag_min_min_sweep, gkcag_min_max_sweep, n_gkcag_steps)
    gkcag_max_range = np.linspace(gkcag_max_min_sweep, gkcag_max_max_sweep, n_gkcag_steps)
    n_gcal_steps = 1
    gcal_min_range = np.linspace(gcal_min_min_sweep, gcal_min_max_sweep, n_gcal_steps)
    gcal_max_range = np.linspace(gcal_max_min_sweep, gcal_max_max_sweep, n_gcal_steps)
    n_gleak_steps = 1
    gleak_min_range = np.linspace(gleak_min_min_sweep, gleak_min_max_sweep, n_gleak_steps)
    gleak_max_range = np.linspace(gleak_max_min_sweep, gleak_max_max_sweep, n_gleak_steps)
    gor_exp = 2

    sim_cnt = 0

    log_file = file('orn_sweep_params.log', 'w')
    for gor_min in gor_min_range:
        for gor_max in gor_max_range:
            for gkcag_min in gkcag_min_range:
                for gkcag_max in gkcag_max_range:
                    for gcal_min in gcal_min_range:
                        for gcal_max in gcal_max_range:
                            for gleak_min in gleak_min_range:
                                for gleak_max in gleak_max_range:
                                    script_cmd = 'python HandTuneOrnParameters.py %d %f %f %f %f %f %f %f %f %f' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max, gor_exp)
                                    log_info = '%d %f %f %f %f %f %f %f %f % f\n' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max, gor_exp)
                                    log_file.write(log_info)
                                    log_file.flush()
#                                    if sim_cnt > 5353:
                                    os.system(script_cmd)
                                    print script_cmd
                                    sim_cnt += 1
                    
                    #gor_min_max = np.array([7e-5, 1e-3])
    #script_cmd = 'python HandTuneOrnParameters.py %d %f %f %f %f %f %f %f %f' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max)


def run_new_parameters(fn=None):

    if fn == None:
        fn = 'new_orn_params.dat'
    d = np.loadtxt(fn)

    n_sim = d[:, 0].size
    log_file = file('orn_sweep_params.log', 'w')
    for i_ in xrange(n_sim):
        sim_cnt = i_
        gor_min = d[i_, 1]
        gor_max = d[i_, 2]
        gkcag_min = d[i_, 3]
        gkcag_max = d[i_, 4]
        gcal_min = d[i_, 5]
        gcal_max = d[i_, 6]
        gleak_min = d[i_, 7]
        gleak_max = d[i_, 8]
        script_cmd = 'python HandTuneOrnParameters.py %d %f %f %f %f %f %f %f %f' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max)
        print 'Command:\n', script_cmd
        log_info = '%d %f %f %f %f %f %f %f %f\n' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max)
        log_file.write(log_info)
        log_file.flush()
        os.system(script_cmd)

if __name__ == '__main__':

#    run_new_parameters(fn='new_combined_orn_params.dat')
#    run_new_parameters(fn='new_orn_params.dat')
    run_full_sweep()

