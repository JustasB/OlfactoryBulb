# -*- coding: utf-8 -*-
import sys
import pylab
import numpy
import math
#from matplotlib import pyplot as plt 
from matplotlib import pylab
from matplotlib import rc # To possibly render in LaTeX
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import lines
import os

class RawLine:
    def __init__(self,src,tgt,sigidx,w):
        self.src=src
        self.tgt=tgt
        self.sigidx=sigidx
        self.w=w

class LoadError(Exception):
    def __init__(self, value='The data has not yet been loaded. Load the data file first.'):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SynWeightSnapshot:
    def __init__(self, nummit=5, numgran=100):
        """Number of dendrodendritic synapses"""
        self.nsyn=0
        
        """Number of mitral cells"""
        self.nummit = nummit
        
        """Number of granule cells"""
        self.numgran = numgran
                
        """Timestamps of each snapshot"""
        self.timestamps=[]
        
        """The mitral-to-granule synaptic weights. Weights are specified in
        3D with rows being the number of mitral cells, columns being the number
        of granule cells, and the z axis is time (that is, a snapshot)."""
        self.m2g=[]
        
        """The granule-to-mitral synaptic weights. Also a 3D structure."""
        self.g2m=[]
        
        """Input data filename, without extension. (Assumes .dat extension)"""
        self.input_filename=None
        
        """Output graphic filename, without the extension."""
        self.output_filename=None
        
        self.data_path=None
        
        """File format"""
        self.file_format='pdf'
        
        """If the output file_format is 'pdf', then this object encapsulates
        multiple pages.
        """
        self.pdf=None

        self.fi_gmax = 0.003

        self.ampa_gmax = 2.
        
    def read_file(self,input_filename,data_path=None,output_filename=None):
        """Read the weight snapshot file and return as a formatted list of line
        objects to be parsed."""
        try:
            fname=''
            if data_path is not None:
                fname=data_path
            if not input_filename.endswith('.dat'):
                input_filename += '.dat'
            filepath = os.path.join(fname, input_filename)
            f=open(filepath,'r') # Open the file
            # The format of the file is such that the first line contains the
            # number of dendrodendritic synapses. The second line contains the
            # timestamp of the first snapshot. 
            # The subsequent lines are formatted as:
            # gid   ->   gid   location_on_sigmoid   value
            # Mitral -> Granule syns will be [0 nummit] in the 1st column, and
            # Granule -> Mitral syns will be [0 nummit] in the second column,
            # with the exception that there will be a timestamp as well that
            # will occur every interval.
            
            raw=[]
            for line in f:
                tosplit=line.split()
                try:
                    c1=int(tosplit[0])
                    c2=int(tosplit[1])
                    c3=int(tosplit[2])
                    c4=float(tosplit[3])
                    raw.append(RawLine(c1,c2,c3,c4))
                except IndexError:
                    # This will occur at the beginning of each snapshot file
                    if self.nsyn==0:
                        self.nsyn=c1
                    else:
                        self.timestamps.append(c1)
                        
                except ValueError:
                    # This will occur if the beginning of the snapshot has a
                    # floating point number as its timestamp.
                    c1=int(math.floor(float(tosplit[0])))
                    self.timestamps.append(c1)

            f.close()
        
        except:
            raise
        
        self.data_path=data_path
        self.input_filename=input_filename[:-4] # Remove '.dat'
        if output_filename is None:
            self.output_filename=self.input_filename+'_img'
        else:
            self.output_filename=output_filename
        return raw

    def parse_data(self,raw):
        """Build the mitral-to-granule 3D volume and the granule to mitral
        3D volume."""
        depth=int(float(len(raw))/float(self.nsyn))
        slice_size=self.nummit*self.numgran*2
        self.m2g=numpy.ndarray(shape=(self.nummit,self.numgran,depth),dtype=float)
        self.g2m=numpy.ndarray(shape=(self.nummit,self.numgran,depth),dtype=float)
        d=0
        count=0
        for r in raw:
            if r.src < self.nummit:
                self.m2g[r.src,r.tgt-self.nummit,d]=r.w
            else:
                self.g2m[r.tgt,r.src-self.nummit,d]=r.w
            count=count+1
            if count % slice_size == 0:
                d=d+1
    
    def get_full_output_path(self,slice=-1):
        fname=''
        if self.data_path is not None:
            fname=self.data_path
        
        if slice==-1 or self.file_format=='pdf':
            fname = os.path.join(fname, '%(a)s.%(c)s' % \
                    {'a':self.output_filename, 'c':self.file_format})
        else:
            fname = os.path.join(fname, '%(a)s%(b)04d.%(c)s' % \
                    {'a':self.output_filename, 'b':slice, 'c':self.file_format})
        return fname
        
    
    def plot_slice(self,slice=0,draw_colorbar=False,space_colorbar=False,cmap='OrRd',figsize=(6,1.5),dpi=144,background_color='gray'):
        # The snapshot slice number to draw, zero-indexed.
        pylab.clf()
        fig=pylab.figure(figsize=figsize) # Make a main figure object
        fig_size=fig.get_size_inches()
        aspect_ratio=float(fig_size[0])/float(fig_size[1])
        # Draw the image such that the excitatory synapses are on top
        # and the inhibitory synapses are on bottom. To do this, we add
        # self.nummit subplots, each subplot containing a subimage
        # containing the excitatory synapse in the first row and the
        # inhibitory synapses in the bottom row.
        m_size=1000 # Won't draw outside, so make it big
        subimg=numpy.zeros((2,self.numgran))
        maxm2g=numpy.max(self.m2g)
        maxg2m=numpy.max(self.g2m)
        ppp=55.    # Short for points per pixel
        tick_label_size=5+int(0.5+sum(fig_size)/aspect_ratio)
        label_font_size=7+int(1.2*sum(fig_size)/aspect_ratio)
        if draw_colorbar is True:
            space_colorbar=True
            
        for i in range(self.nummit):
            ax=fig.add_subplot(self.nummit,1,(i+1)) # Add a new subplot
            subimg[0,:]=self.m2g[i,:,slice]/maxm2g   # Get the mitral-to-granule syn weights
            subimg[1,:]=self.g2m[i,:,slice]/maxg2m   # Get the granule-to-mitral syn weights
            imgplot = ax.imshow(subimg)            # Make an image from the data
            imgplot.set_interpolation('nearest')   # Pixels are not antialiased
            imgplot.set_cmap(cmap)                # Jet colormap
            imgplot.axes.set_aspect('auto')             # Make the pixels taller than wide
            
            # Draw lines between the pixels
            X=range(-1,self.numgran)
            X=numpy.add(X,.5)
            Y=numpy.ones_like(X)*.5
            line=lines.Line2D(X,Y,color=background_color,marker='|',markeredgewidth=1.25,markersize=m_size)
            ax.add_line(line)

            if i<self.nummit-1:                   # Only draw x tick labels on the bottom row
                ax.set_xticklabels([])
            else:
                # Make x labels aligned with where the mitral cells are.
                xticks=numpy.arange(0,self.numgran,float(self.numgran)/float(self.nummit)/2.)
                xticklabels=[]
                for j in range(len(xticks)):
                    if j%2==0:
                        xticklabels.append('')
                    else:
                        xticklabels.append(str(int(round(xticks[j]))))
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels,size=tick_label_size)
            
            # Make 5 y tick marks such that they span correctly to draw the "E" and "I"
            ax.set_yticks(numpy.arange(-.5,2,.5))
            # When horizontally centering, the bounding box of the tick label does not move.
            # position=(-.005,0) moves the text slightly to the left.
            ytick_labels=ax.set_yticklabels(['','E','','I',''], \
            horizontalalignment='center', \
            size=tick_label_size, \
            position=(-.0005*tick_label_size,0))
            ylabel='M'+str(i+1)
            ax.set_ylabel(ylabel, rotation='horizontal', labelpad=label_font_size*.8, size=label_font_size)

            for line in ax.get_yticklines(): # Turn off the Y tickmarks
                line.set_markersize(0)

            imgplot.figure.canvas.draw()    # Command to draw, otherwise the last subimg
                                            #  will be drawn in place of all of them
        left_bound=float(2.5*label_font_size/ppp/fig_size[0])
        
        # Bottom bound will be the height of the font in figure space. 
        # Calculate with 55 points/inch
        bottom_bound=float(1.1*tick_label_size/ppp/fig_size[1])
        top_bound=1-float(1.1*label_font_size/ppp/fig_size[1])
        right_bound=.99
        cb_width=float(.8*tick_label_size/ppp/fig_size[0])
        cb_pad=5./ppp/fig_size[0]
        if space_colorbar is True:
            right_bound-=(cb_width+cb_pad+cb_width)
        fig.subplots_adjust(top=top_bound,bottom=bottom_bound,left=left_bound,right=right_bound)
        
        if draw_colorbar is True:
            cax=fig.add_axes([right_bound+cb_pad,bottom_bound,cb_width,top_bound-bottom_bound])
            a=numpy.outer(numpy.arange(0,1.01,0.01),numpy.ones(1))
            cax.imshow(a,aspect='auto',cmap=cmap,origin="lower")
            yticks=cax.get_ylim()
            cax.set_yticks(yticks)
            # Set position=(2.4,0) because trying to set the position to
            # "right" does not change the font size.
            cax.set_yticklabels(['0','1'],size=tick_label_size,position=(2.4,0))
            cax.set_xticklabels([])
            for line in cax.get_xticklines(): # Turn off the X tickmarks
                line.set_markersize(0)

        titlestr='$t='+str(self.timestamps[slice])+'$' # Make a LaTeX string
        fig.suptitle(titlestr,size=label_font_size)
        fname=self.get_full_output_path(slice)
        #print fname

        if self.file_format=='pdf':
            if self.pdf==None:
                self.pdf=PdfPages(fname)
            fig.savefig(self.pdf, format='pdf', dpi=dpi)
        else:
            fig.savefig(fname,dpi=dpi) 
        pylab.close()
            
    def plot_slices(self,r=[],file_format='pdf',draw_colorbar=False,space_colorbar=False,cmap='OrRd',figsize=(6,1.5),dpi=144,background_color='gray'):
        if len(self.m2g)==0:
            raise LoadError()
        self.file_format=file_format
        self.pdf=None # Let plot_slice() set, if necessary.
        if type(r) is not list:
            r=range(r,r+1) # Assume int, but sage may make it's own integer type.
        if len(r)==0:
            r=range(self.m2g.shape[-1]) # Last parameter of the m2g shape array is the depth
        for i in r:
            self.plot_slice(i,draw_colorbar,space_colorbar,cmap,figsize,dpi,background_color)
            if i%10==0:
                print 'slice '+str(i)
        if self.file_format=='pdf':
            if self.pdf is not None:
                self.pdf.close()
        
        # If using SAGE, the files were probably saved to the DATA folder. 
        # Provide a link in the worksheet to the first and last snapshots for
        # visual output.
        if self.data_path is not None:
            fname=self.get_full_output_path(r[-1])            
            #os.symlink(fname,self.output_filename+'.'+self.file_format)
            #refstr='<a href=\"%(a)s\">%(b)s</a>' % \
            #{'a':fname, 'b':self.output_filename+'.'+self.file_format}
            #html(refstr)
        
sws=SynWeightSnapshot()

def load_file(input_filename,data_path=None,output_filename=None):
    global sws
    raw=sws.read_file(input_filename,data_path,output_filename)
    sws.parse_data(raw)
    
def plot_slices(r=[],file_format='pdf',draw_colorbar=False,space_colorbar=False,cmap='OrRd',figsize=(6,1.5),dpi=144,background_color='gray'):
    global sws
    sws.plot_slices(r,file_format,draw_colorbar,space_colorbar,cmap,figsize,dpi,background_color)
