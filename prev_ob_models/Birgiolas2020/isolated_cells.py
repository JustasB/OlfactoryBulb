from prev_ob_models.utils import RunInClassDirectory, IsolatedCell
import os, sys
import importlib

class OlfactoryBulbCell(IsolatedCell):
    def __init__(self, cell_id):
        self.cell_id = str(cell_id)

        with RunInClassDirectory(OlfactoryBulbCell):
            from neuron import h, load_mechanisms

            # Load the channels - from Mechanisms subfolder
            load_mechanisms("Mechanisms")

            h.load_file("stdrun.hoc")
            h.celsius = 35
            h.cvode_active(1)

            cell_name = self.cell_type + str(cell_id)

            # Load the cell HOC file (they follow MC1.hoc, GC1.hoc, ... pattern)
            os.chdir("Cells")
            self.hoc_path = os.path.abspath(cell_name + ".hoc")
            h.load_file(self.hoc_path)
            os.chdir("..")

            # Build the cell
            self.hoc_template = self.cell_type + str(cell_id)
            self.cell = getattr(h, self.hoc_template)()

            # Apply 3D transformations, if any
            os.chdir("Cells"); sys.path.append(os.getcwd());

            transformation_file = self.cell_type + "Transforms.py"

            if os.path.exists(transformation_file):
                try:
                    # print('Applying: ' + transformation_file)
                    Transform = getattr(importlib.import_module(self.cell_type + 'Transforms'), 'Transform' + cell_name)
                    Transform.apply_on(str(self.cell))
                except:
                    import traceback
                    traceback.print_exc()
                    print('Could not load Transform for '+ cell_name)

            os.chdir("..")

            self.h = h
            self.soma = self.cell.soma

            h.init() # without this, h.run() with multiple cells produces a convergence error

    def set_model_params(self, param_values):
        from neuron import h

        # The diam param functions as a scaling factor of the existing diameters
        # to be able to set the params multiple times, the original diams need to be used when multiplying
        # Here, this is done by instantiating an unmodified copy of the cell
        unmodified_cell = getattr(h, self.hoc_template)()

        for pi, pv in enumerate(param_values):
            attr = self.params[pi]["attr"]
            if attr == "tau_CaPool":
                setattr(h, attr, pv)
            else:
                for param_list in self.params[pi]["lists"]:
                    if not hasattr(self.cell, param_list):
                        # print("List",param_list, "not found in cell", self.hoc_template, "Skipping list...")
                        continue

                    if attr == "diam":
                        unmodified_secs = [sec for sec in getattr(unmodified_cell, param_list)]

                    for isec, sec in enumerate(getattr(self.cell, param_list)):
                        if attr == "diam":
                            for i3d in range(int(h.n3d(sec=sec))):
                                orig_diam = h.diam3d(i3d, sec=unmodified_secs[isec])
                                h.pt3dchange(i3d, orig_diam * pv, sec=sec)
                        elif attr == "L" and param_list == "somatic":
                            L = pv
                            h.pt3dchange(0, -L / 2, 0, 0, L, sec=sec)
                            h.pt3dchange(1,      0, 0, 0, L, sec=sec)
                            h.pt3dchange(2, +L / 2, 0, 0, L, sec=sec)
                        else:
                            setattr(sec, attr, pv)

        del unmodified_cell

class MC(OlfactoryBulbCell):

    # Parameters and their ranges used for fitting
    params = [
        {"attr": "diam", "low": 0.1, "high": 5.0, "lists": ["apical", "basal", "axonal"]},
        {"attr": "L", "low": 6.58, "high": 13.37, "lists": ["somatic"]},

        {"attr": "Ra", "low": 1.0, "high": 150.0, "lists": ["all"]},
        {"attr": "cm", "low": 0.1, "high": 2.0, "lists": ["all"]},
        {"attr": "ena", "low": 20.0, "high": 80.0, "lists": ["all"]},
        {"attr": "ek", "low": -100.0, "high": -50.0, "lists": ["all"]},
        {"attr": "e_pas", "low": -90.0, "high": -50.0, "lists": ["all"]},
        {"attr": "g_pas", "low": 0, "high": 0.00020, "lists": ["all"]},
        {"attr": "sh_Na", "low": 0, "high": 10, "lists": ["all"]},
        {"attr": "tau_CaPool", "low": 1, "high": 300, "lists": ["all"]},

        {"attr": "gbar_Na", "low": 0, "high": 0.2, "lists": ["all"]},
        {"attr": "gbar_Kd", "low": 0, "high": 0.1, "lists": ["all"]},
        {"attr": "gbar_Kslow", "low": 0, "high": 0.002, "lists": ["all"]},
        {"attr": "gbar_KA", "low": 0, "high": 0.02, "lists": ["all"]},
        {"attr": "gbar_KCa", "low": 0, "high": 0.016, "lists": ["all"]},
        {"attr": "gbar_LCa", "low": 0, "high": 0.0005, "lists": ["all"]},

        {"attr": "eh", "low": -40.0, "high": -10.0, "lists": ["apical"]},
        {"attr": "gbar_Ih", "low": 0, "high": 0.000006, "lists": ["apical"]},
        {"attr": "gbar_CaT", "low": 0, "high": 20e-3, "lists": ["apical"]},
    ]

    cell_type = "MC"


class MC1(MC):
    def __init__(self):
        super(MC1, self).__init__(cell_id=1)
        self.param_values = [1.582842487851817, 12.75478806348865, 109.79690032633621, 1.1657185098883014, 22.45075303312652, -61.43158910747635, -69.83937614016634, 3.880041878088015e-05, 6.4220546911647105, 97.95855435014933, 0.09217457722022222, 0.05612925376113354, 0.0008311242905068193, 0.005933006625189661, 5.742606879835478e-05, 0.00038023615700368376, -25.40850883667647, 1.3687634587693794e-06, 0.007113876292779865]
        self.set_model_params(self.param_values)

class MC2(MC):
    def __init__(self):
        super(MC2, self).__init__(cell_id=2)
        self.param_values = [1.0943772088036594, 10.440965950338482, 93.15036142202271, 1.6945282007616136, 25.20516392158557, -74.8777637388441, -59.39666327486791, 3.596078414653267e-05, 2.5554504127540776, 1.389489941154919, 0.0881816004457904, 0.04872161548641748, 0.0003578758191393118, 0.010541253862535379, 0.00847696390171203, 0.0003925972178029332, -19.892718874189566, 7.150418395130422e-07, 0.010948223739908794]
        self.set_model_params(self.param_values)

class MC3(MC):
    def __init__(self):
        super(MC3, self).__init__(cell_id=3)
        self.param_values =  [1.9527098132845895, 11.101410651857844, 22.864583232036104, 0.6427904751592232, 25.205708912391266, -62.37629959504966, -76.93037273748809, 3.174482890928424e-05, 0.010743054467625912, 80.6793136041459, 0.06596349995882334, 0.030788506964382054, 3.938385746162231e-05, 0.0017269117207132935, 0.008743882299355228, 0.0001878679485059271, -33.593831405825284, 5.158712768658259e-07, 0.004634102194662834]
        self.set_model_params(self.param_values)

class MC4(MC):
    def __init__(self):
        super(MC4, self).__init__(cell_id=4)
        self.param_values = [1.1889955393643883, 6.713316586000511, 113.0647154118838, 1.0195884959924766, 27.047779923908816, -70.21445104015474, -58.44571000276015, 2.41179589261979e-05, 2.1234251230802768, 116.66617044665243, 0.08938182311564233, 0.028305145075594434, 0.0013717137679573465, 0.008235801667492574, 0.006295591114313158, 6.059225606189184e-05, -16.431775954888618, 2.5545247291881558e-06, 0.013625818845992767]
        self.set_model_params(self.param_values)

class MC5(MC):
    def __init__(self):
        super(MC5, self).__init__(cell_id=5)
        self.param_values = [1.4659052116669653, 11.272989236805506, 45.81580432923558, 0.6423077604986799, 24.189605615199625, -61.22263984742484, -77.76500118588712, 2.799474934847125e-05, 7.473006638962397, 21.769108516252658, 0.08030263281497557, 0.045378073990493345, 1.940515013623855e-05, 0.003367433087058607, 0.008749055435310703, 0.00011617477848639343, -32.61943576903494, 5.730004860036104e-07, 0.018106878080415785]
        self.set_model_params(self.param_values)


class GC(OlfactoryBulbCell):

    # Parameters and their ranges used for fitting
    params = [
        { "attr": "diam",  "low": 0.1, "high": 3.0, "lists": ["apical"]},
        { "attr": "L",     "low": 0.89, "high": 5.04, "lists": ["somatic"]},

        { "attr": "Ra",    "low": 5.0, "high": 150.0, "lists": ["all"]},
        { "attr": "cm",    "low": 0.1, "high": 10.0, "lists": ["all"]},
        { "attr": "ena",   "low": 10.0, "high": 90.0, "lists": ["all"]},
        { "attr": "ek",    "low": -100.0, "high": -30.0, "lists": ["all"]},
        { "attr": "e_pas", "low": -100.0, "high": -50.0, "lists": ["all"]},
        { "attr": "g_pas", "low": 0, "high": 0.0040, "lists": ["all"]},
        { "attr": "sh_Na", "low": 0, "high": +10, "lists": ["all"]},

        { "attr": "gbar_Na", "low": 0, "high": 0.4, "lists": ["apical"]},
        { "attr": "gbar_Kd", "low": 0, "high": 1.6, "lists": ["apical"]},

        { "attr": "gbar_Na", "low": 0, "high": 5.0, "lists": ["somatic"]},
        { "attr": "gbar_Kd", "low": 0, "high": 5.0, "lists": ["somatic"]},
        { "attr": "gbar_KA", "low": 0, "high": 0.8, "lists": ["somatic"]},
        { "attr": "eh",      "low": -60.0, "high": -10.0, "lists": ["somatic"]},
        { "attr": "gbar_Ih", "low": 0, "high": 0.000200, "lists": ["somatic"]},
        { "attr": "gbar_KM", "low": 0, "high": 0.13, "lists": ["somatic"]},

    ]

    cell_type = "GC"

class GC1(GC):
    def __init__(self):
        super(GC1, self).__init__(cell_id=1)
        self.param_values = [0.6354331329119002, 1.2938842675019608, 86.60228776188752, 6.097757804039867, 25.2960798724282, -76.70759599060172, -89.3885367109972, 0.00011634241105991428, 2.6678571405236333, 0.19655222427059713, 0.20216712905593315, 0.3960861517225469, 1.3914563069980788, 0.506414168406638, -10.019130794184761, 1.5062876653742342e-08, 0.00036501509920944155]
        self.set_model_params(self.param_values)

class GC2(GC):
    def __init__(self):
        super(GC2, self).__init__(cell_id=2)
        self.param_values = [0.5919240225726485, 0.9377868051004455, 135.48423907699538, 9.868734888180338, 44.457228147075426, -69.849717759686, -87.44367826094535, 0.0002532025410260323, 0.6794102899039003, 0.14280039829386304, 0.6754823006187463, 1.5101328902421631, 0.5050189668908865, 0.5196791927599024, -10.798383152114752, 4.855896020236519e-06, 0.08699171893837351]
        self.set_model_params(self.param_values)

class GC3(GC):
    def __init__(self):
        super(GC3, self).__init__(cell_id=3)
        self.param_values = [0.5207279006127108, 2.7418437634170245, 29.86526990214049, 9.915742710613806, 34.811323621837964, -47.61720974003267, -79.25006914322489, 0.00026346830363808586, 8.994121763868439, 0.2060895658072, 0.056796599167294994, 1.298689426862954, 2.8410280349325325, 0.7995603568770102, -58.427640754751636, 7.935912189256338e-07, 0.012519341898903455]
        self.set_model_params(self.param_values)

class GC4(GC):
    def __init__(self):
        super(GC4, self).__init__(cell_id=4)
        self.param_values = [0.2541071033868456, 0.9018686615121411, 5.840902030011179, 8.982427157672912, 33.26154903312686, -72.49039418565542, -87.00745282633099, 0.0002875798281665264, 1.3482155620981673, 0.2329630644884414, 0.7991879969967408, 1.1557906139672314, 1.079595577699155, 0.7422718950734719, -15.525477880615483, 2.407280315875134e-05, 0.0029172228059790853]
        self.set_model_params(self.param_values)

class GC5(GC):
    def __init__(self):
        super(GC5, self).__init__(cell_id=5)
        self.param_values = [0.8779752631390383, 1.858976049156634, 110.66782930144701, 3.9732889416523793, 27.681984895293596, -58.08391616048308, -91.92491812845267, 0.00015430675359767468, 4.8308294713153055, 0.1805187512240916, 1.1277751441096302, 1.989811028917615, 1.155885940569044, 0.3619221096015237, -45.459550085137906, 3.236141764651116e-05, 0.00020697623727848618]
        self.set_model_params(self.param_values)
        
        
        
class TC(OlfactoryBulbCell):

    # Parameters and their ranges used for fitting
    params = [
        {"attr": "diam", "low": 0.1, "high": 2.0, "lists": ["apical", "basal", "axonal"]},
        {"attr": "L", "low": 3.5, "high": 11.6, "lists": ["somatic"]},

        {"attr": "Ra", "low": 1.0, "high": 150.0, "lists": ["all"]},
        {"attr": "cm", "low": 0.1, "high": 5.0, "lists": ["all"]},
        {"attr": "ena", "low": 20.0, "high": 80.0, "lists": ["all"]},
        {"attr": "ek", "low": -100.0, "high": -50.0, "lists": ["all"]},
        {"attr": "e_pas", "low": -90.0, "high": -50.0, "lists": ["all"]},
        {"attr": "g_pas", "low": 0, "high": 0.00040, "lists": ["all"]},
        {"attr": "sh_Na", "low": 0, "high": 10, "lists": ["all"]},
        {"attr": "tau_CaPool", "low": 1, "high": 300, "lists": ["all"]},

        {"attr": "gbar_Na", "low": 0, "high": 0.1, "lists": ["all"]},
        {"attr": "gbar_Kd", "low": 0, "high": 0.2, "lists": ["all"]},
        {"attr": "gbar_Kslow", "low": 0, "high": 0.002, "lists": ["all"]},
        {"attr": "gbar_KA", "low": 0, "high": 0.02, "lists": ["all"]},
        {"attr": "gbar_KCa", "low": 0, "high": 0.016, "lists": ["all"]},
        {"attr": "gbar_LCa", "low": 0, "high": 0.0010, "lists": ["all"]},

        {"attr": "eh", "low": -40.0, "high": -10.0, "lists": ["apical"]},
        {"attr": "gbar_Ih", "low": 0, "high": 0.000060, "lists": ["apical"]},
        {"attr": "gbar_CaT", "low": 0, "high": 20e-3, "lists": ["apical"]},
    ]

    cell_type = "TC"


class TC1(TC):
    def __init__(self):
        super(TC1, self).__init__(cell_id=1)
        self.param_values = [0.5320796104395679, 9.270165432086312, 104.23041806738314, 0.30530014688033, 28.924461902858308, -68.31744332210725, -56.60436604364108, 0.00029150720192176895, 3.913349715740908, 188.34247946378073, 0.08980258012743354, 0.10730949166868958, 0.0003357954496264246, 0.0049970948213131534, 0.0029412988658005113, 1.798259566363729e-05, -38.91275512857027, 7.77139335091131e-06, 0.005384365036192561]
        self.set_model_params(self.param_values)

class TC2(TC):
    def __init__(self):
        super(TC2, self).__init__(cell_id=2)
        self.param_values = [0.7014424046777784, 9.766253820226583, 145.18035124845895, 0.4360545253143458, 28.231448536057798, -66.25059013579343, -69.5821638002925, 0.0002709350365889211, 0.16267963099015992, 10.035443270197067, 0.08923339434172924, 0.1639222113166549, 0.0014699904095443425, 0.007385363730973165, 0.004801901318984596, 0.000498237152921818, -13.055562872667668, 2.491054004831983e-07, 0.012564154351756677]
        self.set_model_params(self.param_values)

class TC3(TC):
    def __init__(self):
        super(TC3, self).__init__(cell_id=3)
        self.param_values = [1.1150775732429254, 4.4342321443458665, 140.0603567489126, 2.5647876077309393, 36.98887466006239, -64.10144100660276, -72.04410633972584, 0.00018889462798377362, 0.9986730975896996, 265.1951143786057, 0.09497953219449135, 0.17075238612729338, 0.0006995167792363354, 0.01258253721755341, 5.983382586310648e-06, 0.0005467916756341761, -12.23003917041239, 1.464842906610013e-05, 0.019510389835192148]
        self.set_model_params(self.param_values)

class TC4(TC):
    def __init__(self):
        super(TC4, self).__init__(cell_id=4)
        self.param_values = [0.880892886671345, 5.427464927311558, 57.73964108732055, 1.2553640380097877, 48.07430525546685, -69.8309330807106, -59.438072156710426, 4.465709904148039e-05, 0.0153674585604624, 118.9260692659365, 0.03672195458892829, 0.1049422156274652, 0.0005348071482369292, 0.001918455378510126, 0.001540229017858929, 0.00034950147735730446, -38.57314310427289, 3.303294139699703e-05, 0.011898977012691088]
        self.set_model_params(self.param_values)

class TC5(TC):
    def __init__(self):
        super(TC5, self).__init__(cell_id=5)
        self.param_values = [0.45739561744658364, 4.43210198915518, 38.306676426103074, 0.338208858836475, 25.78453453122405, -72.19350074054312, -62.7232566025413, 0.00035375861033552975, 0.024296848024877304, 88.40262670945955, 0.09858081240423741, 0.11389061162000566, 0.00134280190081124, 0.014430506220670415, 0.0005401346693437823, 6.938460042954509e-07, -31.55910335355583, 4.860988596316896e-05, 0.010787608496049355]
        self.set_model_params(self.param_values)