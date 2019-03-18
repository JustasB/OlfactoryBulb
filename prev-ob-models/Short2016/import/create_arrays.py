import itertools
import numpy

# below creates symmetrical arrays


B=range(20,780,20)# /1.0 # 0.0 20.0 ... 620.0
#(start, stop, step)


S=range(20,780,20) #/1.0 # [50 50 ... 50 50]

both=list(itertools.product(B,S))

print both

# create asymmetrical array where S<B
#
# go from B=0 to max_B
# for each B value include only S<B
#
#B=numpy.zeros(0,dtype=float,order='C')
#delta_S = 10 # for S
#delta_B = 60 # for B
#max_B = 300
#big_B = B
#big_S = B
#for new_B in range(delta_B, max_B+delta_B, delta_B):
#  for new_S in range(delta_S,new_B+delta_S,delta_S):
#    big_B=numpy.append(big_B,numpy.array([new_B]))
#    big_S=numpy.append(big_S,numpy.array([new_S]))

# override for special cases Shaina requested:

#big_B=[50, 50, 400, 400]
#big_S=[50, 400, 50, 400]
#big_B=[100, 100, 300, 300]
#big_S=[100, 300, 100, 300]

#both=zip(big_B,big_S) # order they appear in the parameter.hoc template
