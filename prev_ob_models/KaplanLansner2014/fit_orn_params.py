import pylab
import numpy as np
from scipy.optimize import leastsq
import simulation_parameters

def gor_func(x, a1, a2, exp):
    """this function calculates the gor for  a given x """
    return a1 * x**exp + a2

class FitOrnParams(object):

    def __init__(self, params):
        self.params = params

    # ---------- gkcag -------------------- 
    def residuals_gkcag(self, p, y, x):
        """
        x: normal x coordinate
        p: parameters of the function to fit, e.g. a and b in y = a * x + b
        """
        err = y - self.peval_gkcag(x, p)
        return err 
        
    def peval_gkcag(self, x, p):
        return p[0] * np.log10(x) + p[1]
    #    return 10 ** p[0] * (x) ** p[1]

    # ------ gcal ------------- 
    def residuals_gcal(self, p, y, x):
        err = y - self.peval_gcal(x, p)
        return err 
        
    def peval_gcal(self, x, p):
        # gcal(gor) = p[0] * log10(gor) + b
        return p[0] * np.log10(x) + p[1]

    # previous trials:
    # gcal(gor) = 10 ** p[0] * gor ** p[1]
#    return p[0] * np.log10(x) + p[1]
#    return 10 ** p[0] * (x) ** p[1]



    # ------ gleak ------------- 
    def residuals_gleak(self, p, y, x):
        err = y - self.peval_gleak(x, p)
        return err 

    def peval_gleak(self, x, p):
        #  gaussian + 1 / gor
    #    return p[0] * np.exp(p[1] * x ** 2) + p[2] / ((x-p[3]) ** p[4])+ p[5]
    #    return p[0] * np.exp(p[1] * x ** 2) + 1 / (1/p[2] + (x - p[3]) ** p[4])+ p[5]
        # simply gaussian
        return p[0] * np.exp(p[1] * x ** 2) + p[2]

    #    return p[0] * x ** p[1] + p[2] * x + p[3]
    #    return (p[0] - p[1] * np.exp( -p[2] * x)) - p[3] * x + p[4]
    #    return p[0] * x * np.exp(p[1] * x)


    def fit_to_list_of_params(self, parameter_array):
        """
        Keyword arguments:
        parameter_array -- 2-dim array, e.g.
                        #"gna","gk","gkcag","gcal","gleak_orn", "tau_cadec"]
        list1 =         [0.5, 5e-2, 1.3e-3, 1.3e-4,     7.75e-5, 1000]
        list2 =         [0.5, 5e-2, 1.35e-3, 1.35e-4,     7.7e-5, 1000]
        list3 =         [0.5, 5e-2, 2.0e-3, 2.0e-4,     7.675e-5, 1000]
        list4 =         [0.5, 5e-2, 3.25e-3, 3.250e-4,     7.65e-5, 1000]
        list5 =         [0.5, 5e-2, 4.3e-3, 4.3e-4,     7.6e-5, 1000]
        list6 =         [0.5, 5e-2, 5.5e-3, 5.5e-4,     7.5e-5, 1000]
        list7 =         [0.5, 5e-2, 6.5e-3, 6.5e-4,     7.0e-5, 1000]
        list8 =         [0.5, 5e-2, 7.5e-3, 7.5e-4,     6.0e-5, 1000]
        list9 =         [0.5, 5e-2, 8.5e-3, 8.5e-4,     5.0e-5, 1000]
        list10=         [0.5, 5e-2, 9.5e-3, 9.5e-4,     4.0e-5, 1000]
        parameter_array = np.array([list1, list2, list3, list4, list5, list6, list7, list8, list9, list10])
        """



        assert (self.params['rel_orn_mit'] == 1), 'Please set params[rel_orn_mit] = 1 in simulation_parameters'
        gor_array = np.zeros(self.params["n_gor"])
        x_min = 0
        x_max = self.params["n_orn_x"]
        gor_exp = self.params["gor_exp"]
        b = self.params["gor_min"]
        a = (self.params["gor_max"] - self.params["gor_min"]) / ((x_max-1)**self.params["gor_exp"] - x_min**self.params["gor_exp"])
        for i in xrange(self.params["n_orn_y"]):  # y-axis
            for j in xrange(self.params["n_orn_x"]):
                orn_id = i * self.params["n_orn_x"] + j
                gor_value = gor_func(j, a, b, gor_exp)
                gor_array[j] = gor_value

        gkcag_list = parameter_array[:,2].tolist()
        gcal_list = parameter_array[:,3].tolist()
        gleak_list = parameter_array[:,4].tolist()

        gor_list = gor_array.tolist()

        # find functions which fit to the found dependencies
        # between the parameter and gor
        # gkcag(gor) = a * log(gor) 
        n_gor = 100
        gor_min = gor_list[0]
        gor_max = gor_list[-1]
        gor_values = np.arange(gor_min, gor_max, (gor_max-gor_min) / n_gor) 


        # fit gkcag
        p0_gkcag = [8.5e-3, 9e-2] # initial guess
        result = leastsq(self.residuals_gkcag, p0_gkcag, args=(gkcag_list, gor_list), maxfev=1000)
        #print "gkcag result", result
        opt_params_gkcag = result[0]
        gkcag_opt = self.peval_gkcag(gor_values, opt_params_gkcag)
        info = "["
        for p in xrange(len(opt_params_gkcag)):
            info += "%.8e, " % (opt_params_gkcag[p])
        info += "%.8e]" % opt_params_gkcag[-1]
        print "gkcag_opt", info

        # gcal(gor) = 10 ** p[0] * gor ** p[1]
        p0_gcal = [1e-3, 1.] # initial: 
        result = leastsq(self.residuals_gcal, p0_gcal, args=(gcal_list, gor_list), maxfev=1000)
        opt_params_gcal = result[0]
        gcal_opt = self.peval_gcal(gor_values, opt_params_gcal)

        info = "["
        for p in xrange(len(opt_params_gcal)):
            info += "%.8e, " % (opt_params_gcal[p])
        info += "%.8e]" % opt_params_gcal[-1]
        print "gcal_opt", info

        # fit parameters for gleak vs gor function
        # a * x * exp(b * x)
        # initial guess for gaussian:
        p0_gleak = [-1e3, -1e6, 1e-6]
        #    return p[0] * x + np.exp(p[1] * x) + p[2]
        result = leastsq(self.residuals_gleak, p0_gleak, args=(np.array(gleak_list), np.array(gor_array)), maxfev=10000)
        opt_params_gleak = result[0]
        info = "["
        for p in xrange(len(opt_params_gleak) - 1):
            info += "%.8e, " % (opt_params_gleak[p])
        info += "%.8e]" % opt_params_gleak[-1]
        print "gleak_opt", info
        gleak_opt = self.peval_gleak(gor_values, opt_params_gleak)

        # gleak with normal scales
        pylab.figure()
        pylab.xlabel("gor")
        gleak_header = "gleak(gor) = %.3e * exp(%.3e * gor ** 2) + %.3e" % (opt_params_gleak[0], opt_params_gleak[1], opt_params_gleak[2])
        output_fn = "gleak_vs_gor.png"
        pylab.ylabel("gleak")
        pylab.plot(gor_array, parameter_array[:,4], '--o')
        pylab.plot(gor_values, gleak_opt)
        pylab.xlim((pylab.xlim()[0] * 0.9, pylab.xlim()[1] * 1.1))
        pylab.ylim((pylab.ylim()[0] * 0.95, pylab.ylim()[1] * 1.05))
        pylab.title(gleak_header)
        yticks = pylab.yticks()[0]
        yticks_strings = ["%.1e" % i for i in yticks]
        pylab.yticks(yticks, pylab.yticks()[0])
        pylab.savefig(output_fn)


        # gleak with log x scale
        output_fn = "gleak_vs_gor_semilogx.png"
        pylab.figure()
        pylab.xlabel("gor")
        pylab.ylabel("gleak")
        gleak_header = "gleak(gor) = %.3e * exp(%.3e * gor^2) + %.3e" % (opt_params_gleak[0], opt_params_gleak[1], opt_params_gleak[2])
        pylab.semilogx(gor_array, parameter_array[:,4], '--o')
        pylab.semilogx(gor_values, gleak_opt)
        pylab.title(gleak_header)
        yticks = pylab.yticks()[0]
        yticks_strings = ["%.1e" % i for i in yticks]
        pylab.yticks(yticks, pylab.yticks()[0])
        pylab.savefig(output_fn)



        pylab.figure()
        pylab.xlabel("gor")
        gkcag_header = "gkcag(gor) = %.3e * log10(gor) + %.3e" % (opt_params_gkcag[0], opt_params_gkcag[1])
        #gkcag_header = "gkcag(gor) = 10^(%.3e) * gor^(%.3e)" % (opt_params_gkcag[0], opt_params_gkcag[1])
        pylab.title(gkcag_header)
        output_fn = "gkcag_vs_gor.png"
        pylab.ylabel("gkcag")
        pylab.plot(gor_array, parameter_array[:,2], '--o')
        pylab.plot(gor_values, gkcag_opt)

        #pylab.loglog(gor_array, parameter_array[:,2], '-o')
        #pylab.loglog(gor_values, gkcag_opt)
        #pylab.semilogx(gor_array, parameter_array[:,2], '-o')
        #pylab.semilogx(gor_values, gkcag_opt)
        pylab.savefig(output_fn)

        pylab.figure()
        pylab.xlabel("gor")
        #gcal_header = "gcal(gor) = 10^(%.3e) * gor^(%.3e)" % (opt_params_gcal[0], opt_params_gcal[1])
        gcal_header = "gcal(gor) = %.3e * log10(gor) + %.3e" % (opt_params_gcal[0], opt_params_gcal[1])
        pylab.title(gcal_header)
        output_fn = "gcal_vs_gor_log10.png"
        pylab.ylabel("gcal")
        #pylab.plot(gor_array, parameter_array[:,3], '-o')
        #pylab.plot(gor_values, gcal_opt)
        pylab.semilogx(gor_array, parameter_array[:,3], '--o')
        pylab.semilogx(gor_values, gcal_opt)
        #pylab.loglog(gor_array, parameter_array[:,3], '-o')
        #pylab.loglog(gor_values, gcal_opt)
        pylab.savefig(output_fn)
        #pylab.ylim((pylab.ylim()[0], pylab.ylim()[1] * 1.1))
        pylab.savefig(output_fn)


    def plot_intedependencies(self, params_gkcag, params_gcal, params_gleak):


        # ---- set the gor values in the most complicated manner ------ 
        gor_array = np.zeros(self.params["n_gor"])
        x_min = 0
        x_max = self.params["n_orn_x"]
        gor_exp = self.params["gor_exp"]
        b = self.params["gor_min"]
        a = (self.params["gor_max"] - self.params["gor_min"]) / ((x_max-1)**self.params["gor_exp"] - x_min**self.params["gor_exp"])
        for i in xrange(self.params["n_orn_y"]):  # y-axis
            for j in xrange(self.params["n_orn_x"]):
                orn_id = i * self.params["n_orn_x"] + j
                gor_value = gor_func(j, a, b, gor_exp)
                gor_array[j] = gor_value
        
        gor_list = gor_array.tolist()
        # find functions which fit to the found dependencies
        # between the parameter and gor
        # gkcag(gor) = a * log(gor) 
        n_gor = 100
        gor_min = gor_list[0]
        gor_max = gor_list[-1]
        gor_values = np.arange(gor_min, gor_max, (gor_max-gor_min) / n_gor) 
        
        gkcag_opt = self.peval_gkcag(gor_values, params_gkcag)
        gcal_opt = self.peval_gcal(gor_values, params_gcal)
        gleak_opt = self.peval_gleak(gor_values, params_gleak)

        fig = pylab.figure()
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        ax1.plot(gor_values, gkcag_opt)
        ax2.plot(gor_values, gcal_opt)
        ax3.plot(gor_values, gleak_opt)

        ax1.set_xlabel('gor')
        ax2.set_xlabel('gor')
        ax3.set_xlabel('gor')
        ax1.set_ylabel('gkcag')
        ax2.set_ylabel('gcal')
        ax3.set_ylabel('gleak')




if __name__ == '__main__':
    param_tool = simulation_parameters.parameter_storage()
    params = param_tool.params
    FOP = FitOrnParams(params)


    list1 =         [0.5, 5e-2, 1.3e-3, 1.3e-4,     7.5e-5, 1000]
    list2 =         [0.5, 5e-2, 1.35e-3, 1.35e-4,   7.4e-5, 1000]
    list3 =         [0.5, 5e-2, 2.0e-3, 2.0e-4,     7.3e-5, 1000]
    list4 =         [0.5, 5e-2, 3.25e-3, 3.250e-4,  7.2e-5, 1000]
    list5 =         [0.5, 5e-2, 4.3e-3, 4.3e-4,     7.1e-5, 1000]
    list6 =         [0.5, 5e-2, 5.5e-3, 5.5e-4,     7.0e-5, 1000]
    list7 =         [0.5, 5e-2, 6.5e-3, 6.5e-4,     7.0e-5, 1000]
    list8 =         [0.5, 5e-2, 7.5e-3, 7.5e-4,     6.0e-5, 1000]
    list9 =         [0.5, 5e-2, 8.5e-3, 8.5e-4,     5.0e-5, 1000]
    list10=         [0.5, 5e-2, 9.5e-3, 9.5e-4,     4.0e-5, 1000]
    param_list = np.array([list1, list2, list3, list4, list5, list6, list7, list8, list9, list10])
    FOP.fit_to_list_of_params(param_list)


    gcal_params =  [4.99086531e-04, 2.26738160e-03, 2.26738160e-03]
    gleak_params = [4.0e-05, -5.18713818e+05, 3.47077557e-05]
#    gkcag_params = [4.99086530e-03, 2.26738160e-02, 2.26738160e-02]
    gkcag_params = [4.99086530e-03, 2.26738160e-02]
    FOP.plot_intedependencies(gkcag_params, gcal_params, gleak_params)
    pylab.show()
