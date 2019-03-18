
# make sym
import custom_params
custom_params.filename = 'g37e1i002'
import params
params.odor_sequence = [ ('Mint', 50, 7050, 6.5e-3) ]
import runsim
runsim.build_part_model([37],[],'g37.dic')
runsim.run()

