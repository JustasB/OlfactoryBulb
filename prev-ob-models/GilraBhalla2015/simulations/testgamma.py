from pylab import *

## delay - table 1 of Carey et al 2009:
## 154 +- 59 ms from 'inspiration start' to 'max of slope times curvature of response'
delay_mean = 154.0e-3
delay_sd = 59.0e-3
## rise-time - table 2 (anaethetized case): 122 +- 32 ms for 10% to 90% rise time
## Cannot use Gaussian for this - it has a long tail!
## but it is correlated with delay; and the easy formula
## for correlated random variables is for Gaussian random variables.
risetime_mean = 121e-3
risetime_sd = 32e-3
delay_risetime_correlation = 0.3 # delay and rise time are positively correlated.
## duration - table 2 (anaesthetized case): 443 +- 119 ms for response above 50% of peak
duration_mean = 443e-3
duration_sd = 119e-3

def gamma_scale_shape(mean,sd):
    theta = sd**2/mean
    return mean/theta, theta

scale,shape = gamma_scale_shape(delay_mean,delay_sd)
print scale,shape
figure()
gammas = gamma(scale,shape,100000)
hist(gammas,100)
title('delay distrib')
scale,shape = gamma_scale_shape(risetime_mean,risetime_sd)
print scale,shape
figure()
gammas = gamma(scale,shape,100000)
hist(gammas,100)
title('risetime distrib')
scale,shape = gamma_scale_shape(duration_mean,duration_sd)
print scale,shape
figure()
gammas = gamma(scale,shape,100000)
hist(gammas,100)
title('duration distrib')


k,theta = gamma_scale_shape(4.0,sqrt(8.0))
print k,theta
figure()
gammas = gamma(k,theta,100000)
hist(gammas,100)
show()
