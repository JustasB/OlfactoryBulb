
# coding: utf-8

# In[1]:


import os; os.chdir("..");
import multiprocessing
from prev_ob_models.Birgiolas2020.fitting import *


# In[3]:


import sys
cell_id = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1].isdigit() else 4
print(cell_id)


# In[4]:


fitter = CellFitter(cell_type="mc", fitting_model_class="prev_ob_models.Birgiolas2020.isolated_cells.MC"+str(cell_id))


# In[ ]:


#[p for p in fitter.params if p["attr"] == "gbar_KCa"][0]["high"] = 0.008


# In[ ]:


fitter.clear_cache()


# In[ ]:


fitter.fit(20, 30)


# In[ ]:


print(fitter.best)


# In[ ]:


fitter.parameter_report(fitter.best)


# In[ ]:


raise


# In[ ]:


from scipy.optimize import minimize
def evaluate(ind):
    return fitter.evaluate(ind)[0]

min_res = minimize(evaluate, fitter.best, method='Powell', options={"disp":True, "maxfev":1})
min_res


# In[5]:


fitter.clear_cache()
#import neuronunit
#neuronunit.tests.plot = True
#fitter.best = [1.284354439836445, 143.98536456930907, 0.688826897386708, 27.447507623212005, -63.2076648832886, -79.13471472057435, 2.5485826914427995e-05, 4.195388939600915, 5.2947644875144375, 0.02557986838703772, 0.018719998635515288, 0.0008312741862146728, 0.00289725491018908, 0.0037058811189542794, 0.0003180609451649635, -26.395680618374698, 1.8686297649494659e-06, 0.0074009036132169775]
fitter.best = []
df, score = fitter.get_best_score()
print(score["model_score"]) 
df
#raise


# In[ ]:


score


# In[ ]:


fitter.clear_cache()
df = fitter.get_previous_model_scores()
df


# In[ ]:


df.to_csv("ephyz/"+fitter.cell_type+"_model_zscores.csv")


# In[ ]:


raise


# In[ ]:


np.linspace(0,100,0.1)


# In[ ]:


score["model_score"]


# In[ ]:


raise


# In[ ]:


fitter.top = None


# In[ ]:


cpus = multiprocessing.cpu_count()
fitter.fit(int(cpus*1.5), 48)


# In[ ]:


fitter.clear_cache()
fitter.best = []
df, score = fitter.get_best_score()
df


# In[ ]:


score["model_score"]


# In[ ]:


fitter.pretty_pop()


# In[ ]:


from prev_ob_models.Birgiolas2020.isolated_cells import MC1
mc = MC1()


# In[ ]:


fitter.print_params(fitter.best, mc.cell)

