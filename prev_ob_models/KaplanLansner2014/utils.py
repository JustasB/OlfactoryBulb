import numpy as np
import os
import re



def distribute_n(n, n_proc, pid):
    """
    l: list of elements to be distributed among n_proc processors
    pid: (int) process id of the process calling this function
    n_proc: total number of processors
    Returns the min and max index to be assigned to the processor with id pid
    """
    n_per_proc = int(n / n_proc)
    R = n % n_proc
    offset = min(pid, R)
    n_min = int(pid * n_per_proc + offset)
    if (pid < R):
        n_max = int(n_min + n_per_proc + 1)
    else:
        n_max = int(n_min + n_per_proc)
    return (n_min, n_max)


def get_recorded_gids(folder, fn_base, ending='.v'):
    """
    Returns a list of cell gids which have been found as last digit in the filenames
    stored in folder and matching fn_base.
    e.g.
    folder = 'Results/VoltageTraces/'
    fn_base = 'pyr_volt_0_' 
    it returns all pyramidal cell gids that have been recorded during this pattern
    """
    gids_found = []
    for fn in os.listdir(folder):
        m = re.match('%s(\d+)(.*)' % fn_base, fn)
        if m:
            path = folder + '/' + fn
            gid = m.groups()[0]
            gids_found.append(int(gid))
    return gids_found


def sort_gids_by_nspikes(fn, gid_offset=0):
    """
    Returns the cell_gids and the corresponding number of spikes fired by those cells
    sorted by the number of spikes in ascending order.
    """

    nspikes = np.loadtxt(fn)
    sorted_idx = nspikes[:, 1].argsort()
    idx_ = nspikes[sorted_idx, 1].nonzero()[0]
    sorted_gids = np.array(nspikes[sorted_idx[idx_], 0], dtype=np.int)
    sorted_nspikes = np.array(nspikes[sorted_idx[idx_], 1], dtype=np.int)
    return sorted_gids, sorted_nspikes

