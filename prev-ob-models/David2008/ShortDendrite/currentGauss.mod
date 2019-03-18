TITLE Internal and synaptic noise


: Include internal noisy channels with density form
: With cut off frequency

: in=-norm(mean,std)

: written by Francois David and Yi Zhou 2005


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS  current_gauss
	
	NONSPECIFIC_CURRENT  in  : negative current depolarizes the membrane
	RANGE del,dur
        
	RANGE in
	RANGE rand
	RANGE mean, std0, std
	RANGE f0  : the sampling frequency
	RANGE N_smooth : =1000/f0/dt
	RANGE count
	RANGE noise_seed
	RANGE tau_f	
	RANGE indic_max,indic_sim,indic_kern
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	del=0 (ms)
	dur=100 (ms) <0,1e9>
	dt    = 0.01 (ms) 
	mean=0         (nA)  : default
	std0=1e-4       (nA)  : default
	std=1e-4
	noise_seed = 1
    
	f0=4000   (1/s)   : default 4000Hz
	tau_f = 5 (ms)    : low pass filter
	t_kern = 20 (ms)  : length of filter
	t_sim = 1000 (ms) : length of stimulus 
}

ASSIGNED {
	v      (mV)

	in	(nA)
	in_temp	(nA)
	rand
	count
        N_smooth
	noise[29000]
	filter[1000]
	current[28000]
	indic_kern
	indic_sim
	indic_max
	uu

	current_total
	noise_total
	current_mean
	noise_mean
	current_var
	noise_var
	power_ratio
}

PROCEDURE seed(x) {
	set_seed(x)
}
 
LOCAL u,w

BREAKPOINT {

	if (t<del+dur-.1 && t>=del) {
		if(count/N_smooth==1){
			in = -current[uu]
			uu=uu+1
			count = 0
		}else{
			count = count+1
                        
		}
	} else {
	  in=0 
		:printf("uu=%g\n",uu)
	}
:printf("count=%g\n",count)
}

INITIAL {

	N_smooth=floor(1000/f0/dt)*2  : break point called twice,here 1000 is a scalar

		:### equalize the power of the sampled white noise2
		:### std_a/std_b=sqrt(f_a/f_b) with f_a=4000 and std0=std_a : see notebook p186

	std=std0/sqrt(4000/f0)  
        	:###note that for Wiener process should be std_a/std_b=sqrt(f_b/f_a)
	
	rand = normrand(mean, std)
	in  = -rand

	count=0
	uu=0
	w =0
	u = 0
 	power_ratio = 0

	indic_kern = t_kern*f0/1000
	indic_sim = dur*f0/1000
	indic_max = indic_sim+indic_kern  : useful for the the convolution product

	FROM u=0 TO indic_kern-1 {
        	filter[u] = exp(-u*(1000/f0)/tau_f)
	}

	FROM u=0 TO indic_max-1 {
		noise[u] = normrand(mean, std)
	}

	FROM u=0 TO indic_sim-1 {
		current[u] = 0
		FROM w=0 TO indic_kern-1 {
			current[u] = current[u] + noise[w+u] * filter[indic_kern-w-1] 
		}
	}

	current_total=0
	noise_total=0

	FROM u=0 TO indic_sim-1 {
		current_total = current_total + current[u]
		noise_total   = noise_total + noise[u]
	}

	current_mean = current_total/(dur*f0/1000)
	noise_mean = noise_total/(dur*f0/1000)
	printf("current_mean= %g\n",current_mean)
	printf("noise_mean= %g\n",noise_mean)


	current_var=0
	noise_var = 0
	FROM u=0 TO indic_sim-1 {
		current_var = current_var + pow(current[u]-current_mean,2)
		noise_var = noise_var + pow(noise[u]-noise_mean,2)
	}
	current_var = current_var/(dur*f0/1000)
	noise_var = noise_var/(dur*f0/1000)

	:printf("current_var1= %g\n",current_var)
	:printf("noise_var1= %g\n",noise_var)
	power_ratio = current_var/noise_var
	:printf("power_ratio= %g\n",power_ratio)

	:Normalization
	FROM u=0 TO indic_sim-1 {
		current[u] = (current[u]-current_mean+noise_mean)/sqrt(current_var/noise_var)
	}
	count = 0

	:OPTIONS -Checking
	:	current_total=0
	:	noise_total=0
	:	FROM u=0 TO indic_sim-1 {
	:		current_total = current_total + current[u]
	:		noise_total   = noise_total + noise[u]
	:	}
	:	current_mean = current_total/(dur*f0/1000)
	:	noise_mean = noise_total/(dur*f0/1000)
	:	printf("current_mean= %g\n",current_mean)
	:	printf("noise_mean= %g\n",noise_mean)
	:	current_var=0
	:	noise_var = 0
	:	FROM u=0 TO indic_sim-1 {
	:		current_var =current_var+ pow(current[u]-current_mean,2)
	:		noise_var = noise_var+pow(noise[u]-noise_mean,2)
	:	}
	:	current_var = current_var/(dur*f0/1000)
	:	noise_var = noise_var/(dur*f0/1000)
	:	printf("current_var2= %g\n",current_var)

	:	printf("noise_var2= %g\n\n",noise_var)
	:	power_ratio = current_var/noise_var
	:	printf("power_ratio2= %g\n",power_ratio)
}

UNITSON
