import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import pylab
import numpy as np
import sys
import os
import MergeSpikefiles 
import simulation_parameters
from ExponentialFiltering import filter_spike_train

import matplotlib
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy.random as rnd
import pylab
import os
import random
from matplotlib.collections import PatchCollection



class CalciumImagingTool(object):
        
    def __init__(self, pn=0):
        self.pn = pn
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)#, axisbg='k')
        rnd.seed(1)
        random.seed(1)
        params2 = {'backend': 'eps',
                  'axes.labelsize': 14,
                  'text.fontsize': 14,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'legend.pad': 0.2,     # empty space around the legend box
                  'legend.fontsize': 14,
                  'lines.markersize': 3,
                  'font.size': 14,
                  'path.simplify': False,
                  'figure.figsize': self.get_figsize(800)}
        pylab.rcParams.update(params2)

    def set_figsize(self, fig_width_pt):
        pylab.rcParams['figure.figsize'] = self.get_figsize(fig_width_pt)

    def get_figsize(self, fig_width_pt):
        inches_per_pt = 1.0/72.0                # Convert pt to inch
        golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        fig_height = fig_width*golden_mean      # height in inches
        fig_size =  [fig_width,fig_height]      # exact figsize
        return fig_size


    def create_cell_paths_regular(self, n_cells):
        scale = 0.3 # shrink factor for cell center placement (cell x, y is at scale * (x,y))
        self.x_max = 30 
        self.y_max = n_cells / self.x_max + 1

        cell_path_dict = {}
        cells = range(n_cells)
        for cell in xrange(n_cells):
            x = cell % self.x_max
            y = int( cell / self.y_max)
#            center = (x + (0.5 - rnd.rand())/3., y + (0.5 - rnd.rand())/3.)
            center = (x + (0.5 - rnd.rand()), y + (0.5 - rnd.rand()))
    #        center = (x + rnd.rand(), y + rnd.rand())
            coords = self.set_npoints(center, 6)
        #    coords = set_points((x, y))
            for i in xrange(int(x)):
                coords = self.rotate_list(coords)
            pathdata = self.build_pathdata(coords)
            codes, verts = zip(*pathdata)
            path = mpath.Path(verts, codes)
            # assign the cell to which these paths belong to
            cell_path_dict[cell] = path 

        return cell_path_dict

    def create_cell_paths(self, n_cells):
        self.x_max = np.sqrt(n_cells)
        self.y_max = self.x_max

        cell_path_dict = {}
        cells = range(n_cells)
        for cell in xrange(n_cells):
            x = cell % self.x_max
            y = int( cell / self.y_max)
            center = (x + (0.5 - rnd.rand())/3., y + (0.5 - rnd.rand())/3.)
    #        center = (x + rnd.rand(), y + rnd.rand())
            coords = self.set_npoints(center, 6)
        #    coords = set_points((x, y))
            for i in xrange(int(x)):
                coords = self.rotate_list(coords)
            pathdata = self.build_pathdata(coords)
            codes, verts = zip(*pathdata)
            path = mpath.Path(verts, codes)
            # assign the cell to which these paths belong to
            rnd_cell = random.choice(cells)
            cells.remove(rnd_cell)
            cell_path_dict[rnd_cell] = path 

        return cell_path_dict


    def set_cells(self, ca_data_fn):
        self.ca_data = np.loadtxt(ca_data_fn)
        # ca_data = np.array((n_cells, n_images + 1))
        # ca_data[:,0] = cell_ids
        # ca_data[gid, i] = Calcium concentration of cell with gid at time i
        
        cells = self.ca_data[:,0]
        n_cells = cells.size
        # draw each cell, i.e. get a dict {gid : path_patch_representing_the_cell}
#        self.cell_path_dict = self.create_cell_paths_regular(n_cells)
        self.cell_path_dict = self.create_cell_paths(n_cells)


    def make_ca_images(self, ca_data_fn, pn=0):    
        self.ca_data = np.loadtxt(ca_data_fn)
        n_cells = self.ca_data[:,0].size
        max_value = np.max(self.ca_data[:,1:])
        min_value = np.min(self.ca_data[:,1:])
        mean_value = self.ca_data[:,1:].mean()
        n_images = self.ca_data[0,1:].size

        baseline = min_value
    #    baseline = self.ca_data[:,1:].mean()
        print "min", baseline
        print "max", max_value
        print "mean", self.ca_data[:,1:].mean()
        patches = [] 
    #    norm = matplotlib.mpl.colors.Normalize(vmin=0, vmax=max_value)
        norm = matplotlib.mpl.colors.Normalize(vmin=0, vmax=2 * mean_value)

        thresh_frac = 0.2
#        thresh = (max_value-min_value) * thresh_frac
        thresh = mean_value
        alpha_mult = 0.82
#        for image in xrange(n_images - 1,n_images):
        for image in xrange(n_images):
            print "Plotting image", image ," for pattern ", pn
            patches = [] 
            values = []
            cells_above_thresh = 0
            for cell in self.cell_path_dict.keys():
                # cell - offset ? 
    #            alpha = self.ca_data[cell, image+1] / max_value

#                alpha = (self.ca_data[cell, image+1] - baseline * 0.9) / (max_value - baseline * 0.9)
                alpha = (self.ca_data[cell, image+1] - baseline) / (max_value - baseline)

#                alpha *= (1 - np.exp(alpha_mult * alpha) + 1)
                value = self.ca_data[cell, image+1] - baseline
                if (self.ca_data[cell, image+1] > thresh):
                    cells_above_thresh += 1
                values.append(value)
    #            print cell, value, self.ca_data[cell, image]
                    
                # draw a patch for the cell with the paths corresponding to the cell
                path = self.cell_path_dict[cell]
    #            patch = mpatches.PathPatch(path, facecolor='red', edgecolor='yellow', alpha=alpha)
#                patch = mpatches.PathPatch(path, facecolor='red', edgecolor="None", alpha=value)
#                if pn == 1:
#                    patch = mpatches.PathPatch(path, facecolor='red', edgecolor="None", alpha=alpha)
#                if pn == 16:
#                    patch = mpatches.PathPatch(path, facecolor='blue', edgecolor="None", alpha=alpha)
                patch = mpatches.PathPatch(path, facecolor='red', edgecolor="None", alpha=alpha)
    #            patch = mpatches.PathPatch(path)#, alpha=alpha)
                patches.append(patch)
                self.ax.add_patch(patch)

            # save the image
#            self.ax.set_xticks(())
#            self.ax.set_yticks(())
            self.ax.set_xlim((-1, self.x_max+1))
            self.ax.set_ylim((-1, self.y_max+1))
#            self.ax.set_ylim((-1, (self.y_max+1) * 0.3))
            self.ax.set_title('Artificial calcium image of \ncortical activity in response to odorant %d' % pn)
            print "In image %d %d / %d  (%f percent) cells are above the activity threshold of %f (%.2f percent)" % (image, cells_above_thresh, n_cells, float(cells_above_thresh) / n_cells * 100, thresh, thresh_frac * 100)
            
    #        plt.show()
                
            # optional colormap
    #        collection = PatchCollection(patches, cmap=matplotlib.cm.jet, norm=norm, alpha=0.6)
    #        collection = PatchCollection(patches, cmap=matplotlib.cm.jet, norm=norm, alpha=0.6)
            collection = PatchCollection(patches)
            collection.set_array(np.array(values))
    #        self.ax.add_collection(collection)
    #        pylab.colorbar(collection)

            fn = "ca_image_%02d_%02d.png" % (pn, image)
            pylab.savefig(fn)
#            fn = "ca_image_%d_%d.svg" % (pn, image)
#            pylab.savefig(fn)
#            self.fig.clear()

    def rotate_list(self, l):
        l2 = list(l)
        last = l2.pop()
        l2.reverse()
        l2.append(last)
        l2.reverse()
        return l2

    def build_pathcoords(self, l):
        # l = list
        # put the first element on the second last position
        # return list
        l2 = list(l)
        last = l2.pop()
        l2.append(l2[0])
        l2.append(last)
        return l2

    def set_npoints(self, center, n_points):
        # n_points 
        r = 0.25  # radius
        assert ((n_points > 2))
        angle_step = 2 * np.pi / n_points
        a = 0   # angle
        coords = []
        for p in xrange(n_points):
            y = center[1] - r * np.cos(a) - r * (rnd.rand() - 0.5) * np.sign(np.sin(a))
            x = center[0] - r * np.sin(a) - r * (rnd.rand() - 0.5) * np.sign(np.sin(a))
            a += angle_step
            coords.append((x,y))
        return coords

    def set_points(self, c):
        # c = center 
        # r = radius
        r = 3.0
        coords = [    (c[0] - r - rnd.rand(),       c[1] - r - rnd.rand())]     # 0
        coords.append((c[0] - 2 * r + rnd.rand(),   c[1] - r/2. + rnd.rand()))  # 1
        coords.append((c[0] - r/2. - rnd.rand(),    c[1] + r + rnd.rand()))     # 2
        coords.append((c[0] + r/2. + rnd.rand(),    c[1] + r + rnd.rand()))     # 3
        coords.append((c[0] + r + rnd.rand(),       c[1] - r/2. + rnd.rand()))  # 4
        coords.append((c[0] + r/2. + rnd.rand(),    c[1] - r - rnd.rand()))     # 5
        return coords


    def build_pathdata(self, l):
        # l = coords (list of (x,y) tuples)
    #    last = l.pop()
    #    l.append(l[0])
    #    l.append(last)
        pathdata = [(Path.MOVETO, l[0])]
        for c in xrange(1, len(l)):
            pathdata.append((Path.CURVE3, l[c]))
        pathdata.append((Path.CURVE3, l[0]))
        return pathdata

    def make_movie(self, png_filter, pn):
        # make the movie of rasterplots
        calcium_movie = "calcium_movie_pattern%d.mp4" % pn
        fps = 2
        os.system("rm %s" % calcium_movie)
        command = "ffmpeg -f image2 -r %f -i %s -b 36000 %s" % (fps, png_filter, calcium_movie)
        os.system(command)
        print calcium_movie


if __name__ == '__main__':

    # MAIN
    n_patterns = 25
    data_path = "/home/bernhard/workspace/Poster_Neurochem/data/calcium_data/"
    #for pn in xrange(3, n_patterns):
    #for pn in xrange(n_patterns):
    for pn in [3, 4, 5]:
        ca_data_fn = data_path + "calcium_mean_data_" + str(pn) + ".dat"
        CIT = CalciumImagingTool()
        CIT.set_cells(ca_data_fn)
        CIT.make_ca_images(ca_data_fn, pn)
        png_filter = "ca_image_%02d" % pn + "_%02d.png"
        print png_filter
        CIT.make_movie(png_filter, pn)
    #    del CIT

