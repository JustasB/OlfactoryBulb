import os
import sys
import time
import simulation_parameters # defines simulation parameters
# classes for setting up connectivity and the individual cell parameters
import CreateObConnections
import CreateOrnParameters


def prepare_epth_ob(params):

    print "Creating ORN parameters...."
    OrnParamClass = CreateOrnParameters.CreateOrnParameters(params)
    ok = OrnParamClass.create_single_odorant_patterns()

    # EPTH -> OB connections are not affected by the pattern
    #print "Creating connections: orn -> mit"
    ConnectionClass = CreateObConnections.CreateObConnections(params)
    ConnectionClass.connect_orn_mit()
    ConnectionClass.connect_orn_pg()
    ConnectionClass.connect_pg_mit_serial()
    ConnectionClass.connect_pg_mit_reciprocal()
    ConnectionClass.connect_mt_gran_local()
    ConnectionClass.connect_mt_gran_global()


if __name__ == '__main__':
    t1 = time.time()

    param_tool = simulation_parameters.parameter_storage()
    # params is the dictionary with all parameters
    params = param_tool.params
    param_tool.write_parameters_to_file(params["info_file"]) # human readable
    param_tool.write_parameters_to_file() # 
    param_tool.hoc_export() # write the simulation parameters to a NEURON executable file
    param_tool.print_cell_gids()

    prepare_epth_ob(params)

    print "Folder name:", params['folder_name']
    t2 = time.time()
    print "Time: %.1f sec %.1f min" % (t2 - t1, (t2-t1)/60.)
