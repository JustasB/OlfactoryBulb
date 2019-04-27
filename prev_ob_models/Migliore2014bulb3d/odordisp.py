from colors import palette as nrncolors
try:
    from enthought.traits.ui.api import Handler
    from enthought.traits.ui.api import UIInfo
except:
    from traitsui.api import Handler
    from traitsui.api import UIInfo
class OdorHandler(Handler):
    
    def __init__(self, fig, bulb, glomval):
        from copy import copy
        self.__glomval = copy(glomval)
        self.__fig = fig; self.__bulb = bulb

        # create method for every odor
        code = "def _show%g(self, info): self.show(%g)"
        for i in range(len(self.__bulb.real_h)):
            x = code % (i, i)
            compcode = {}
            exec x.strip() in compcode
            name = '_show%g' % i
            setattr(self.__class__, name, compcode[name])
            
    # clean the glomerulus
    def clean(self):
        for i in range(len(self.__bulb.real_h)):
            self.__bulb.real_h[i].property.color = (1., 0., 0.)
        self.__fig.scene.render()


    def _clean(self, info): self.clean()
    def _setview(self, info): self.setview()

    # set best view for the camera
    def setview(self):
        # set the camera
        self.__fig.scene.camera.position = [ -294.68175837,  2301.60733878,  5760.60463821 ]
        self.__fig.scene.camera.focal_point = [  972.90044282,  1081.82497137,   -90.96137608 ]
        self.__fig.scene.camera.view_angle = 30.
        self.__fig.scene.camera.view_up = [ 0.04971769,  0.97984323, -0.19348226 ]
        self.__fig.scene.render()
        
    # color gloms
    def show(self, i):
        for j, x in enumerate(self.__glomval[i]):
            self.__bulb.real_h[j].property.color = nrncolors[x]
        self.__fig.scene.render()

def OdorsInput(fname):
    odorlbl = []; glomval = []
    f = open(fname, 'r')

    minval = 100.; maxval = -1.

    line = f.readline()
    while line:
        token = line.split()

        # labels
        lbl = token[0].replace('_', ' ')
        odorlbl.append(lbl)

        # glom values
        vals = []
        for i in range(1, len(token)):
            vals.append(float(token[i]))
        _min = min(vals)
        _max = max(vals)
        if _min < minval: minval = _min
        if _max > maxval: maxval = _max
        glomval.append(vals)
        
        line = f.readline()
    f.close()


    # normalization for a best visualization of colors
    for i in range(len(glomval)):
        for j in range(len(glomval[i])):
            glomval[i][j] -= minval
            glomval[i][j] /= (maxval - minval)
            glomval[i][j] = int(glomval[i][j] * (len(nrncolors) - 1))
    return odorlbl, glomval

def initOdorsDisp(fname, fig, bulb):
    odorlbl, glomval = OdorsInput(fname)
    try:
        from enthought.traits.ui.menu import Action, MenuBar, Menu, Separator # create odor list
    except:
        from traitsui.menu import Action, MenuBar, Menu, Separator # create odor list


    menu = Menu(name='Odors') 
    for i, name in enumerate(odorlbl): menu.append(Action(name=name, action='_show%g' % i))
    #menu.append(Separator())
    menu1 = Menu(name='View') 
    menu1.append(Action(name='Set View as Vinci\'s', action='_setview'))
    menu1.append(Action(name='Clean gloms', action='_clean'))
    return MenuBar(menu, menu1), OdorHandler(fig, bulb, glomval)

