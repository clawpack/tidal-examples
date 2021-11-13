
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from clawpack.geoclaw import topotools
from clawpack.visclaw import gaugetools, plottools, colormaps

cmin_water = -1.
cmax_water = 1.
cmin_speed = 0.
cmax_speed = 2.

def speed(current_data):
    """
    Return a masked array containing the speed only in wet cells.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    from clawpack.visclaw.geoplot import drytol_default
    drytol = getattr(current_data.user, 'drytol', drytol_default)
    q = current_data.q
    h = q[0,:,:]
    hu = q[1,:,:]
    hv = q[2,:,:]
    eta = q[3,:,:]

    hwet = np.ma.masked_where(h <= drytol, h)
    speed = np.sqrt(hu**2 + hv**2) / hwet  # will be masked where dry

    try:
        # Use mask covering coarse regions if it's set:
        m = current_data.mask_coarse
        speed = numpy.ma.masked_where(m, speed)
    except:
        pass

    return speed

if 0:
    sea_level = -0.468
    tide_data = np.loadtxt('tidedata_dec2015.txt', skiprows=1)
    t_tide = tide_data[:,0] / 3600.
    #offset = 2.561  # from MLLW to MHW
    offset = 0.
    eta_tide = interp1d(t_tide, tide_data[:,1]-offset, bounds_error=False,
                        fill_value = np.nan)

if 1:
    # forcing:
    sea_level = 0.
    tperiod = 12  # hours
    eta_tide = lambda t: 1.*np.sin(2*np.pi*t/tperiod)
    t_tide = np.linspace(0,48,1000)



def plot_tide(current_data):
    from pylab import linspace,plot,grid,title,cos,pi,xlabel,ylabel,xticks
    from pylab import xlim,ylim
    tframe = current_data.t / 3600.  
    tt = linspace(t_tide[0],t_tide[-1],1000)
    plot(tt, eta_tide(tt), 'k')
    plot([tframe],[eta_tide(tframe)], 'ro')
    grid(True)
    title('Tide stage in ocean')
    xticks(np.arange(0,5.1*24,24))
    xlabel('hours')
    ylabel('meters (0=MTL)')
    xlim(0,48)
    ylim(-1.2,1.2)

def plot_src_boxes():

    if 0:
        # Johns River
        x1rs = -123.951
        x2rs = -123.950
        y1rs = 46.8815
        y2rs = 46.8825

    if 0:
        # Chehalis River:
        x1rs = -123.743
        x2rs = -123.740
        y1rs = 46.956
        y2rs = 46.958
        kwargs = {'color':'yellow', 'linewidth':1}
        plottools.plotbox([x1rs,x2rs,y1rs,y2rs],kwargs)

        # Wishkah River
        x1rs = -123.811
        x2rs = -123.810
        y1rs = 46.990
        y2rs = 46.991
        plottools.plotbox([x1rs,x2rs,y1rs,y2rs],kwargs)

        # Hoquiam River
        x1rs = -123.886
        x2rs = -123.885
        y1rs = 46.995
        y2rs = 46.996
        plottools.plotbox([x1rs,x2rs,y1rs,y2rs],kwargs)


#--------------------------
def setplot(plotdata=None):
#--------------------------

    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'    # 'ascii' or 'binary' to match setrun.py


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=[0,1102,1187], format_string='ko', add_labels=True)
    

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Domain', figno=50)
    plotfigure.kwargs = {'figsize':(13,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.axescmd = 'axes([.1,.1,.85,.6])'
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = [-124.65,-123.65]
    plotaxes.ylimits = [46.8, 47.1]

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Surface at %4.2f hours' % t, fontsize=15)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        pylab.gca().set_aspect(1/pylab.cos(47*pylab.pi/180))
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin_water
    plotitem.pcolor_cmax = cmax_water
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.colorbar_kwargs = {'extend':'both'}
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    plotaxes = plotfigure.new_plotaxes('1d_plot')
    plotaxes.axescmd = 'axes([.1,.75,.7,.15])'
    plotaxes.afteraxes = plot_tide
    
    

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.kwargs = {'figsize':(13,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.axescmd = 'axes([.1,.1,.85,.6])'
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    #plotaxes.xlimits = [-124.23,-123.65]
    plotaxes.xlimits = [-124.25,-123.65]
    plotaxes.ylimits = [46.82, 47.06]

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Surface at %4.2f hours' % t, fontsize=15)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        pylab.gca().set_aspect(1/pylab.cos(47*pylab.pi/180))
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin_water
    plotitem.pcolor_cmax = cmax_water
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.colorbar_kwargs = {'extend':'both'}
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    plotaxes = plotfigure.new_plotaxes('1d_plot')
    plotaxes.axescmd = 'axes([.1,.75,.7,.15])'
    plotaxes.afteraxes = plot_tide
    
    

    #-----------------------------------------
    # Figure for speed
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Speed', figno=1)
    plotfigure.kwargs = {'figsize':(13,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.axescmd = 'axes([.1,.1,.85,.6])'
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    #plotaxes.xlimits = [-124.23,-123.65]
    plotaxes.xlimits = [-124.25,-123.65]
    plotaxes.ylimits = [46.82, 47.06]

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Speed at %4.2f hours' % t, fontsize=15)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        pylab.gca().set_aspect(1/pylab.cos(47*pylab.pi/180))
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = speed
    plotitem.pcolor_cmap = colormaps.blue_yellow_red
    plotitem.pcolor_cmin = cmin_speed
    plotitem.pcolor_cmax = cmax_speed
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.6
    plotitem.colorbar_kwargs = {'extend':'max'}
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    plotaxes = plotfigure.new_plotaxes('1d_plot')
    plotaxes.axescmd = 'axes([.1,.75,.7,.15])'
    plotaxes.afteraxes = plot_tide

    #-----------------------------------------
    # Figure for zoom
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Aberdeen', figno=10)
    plotfigure.kwargs = {'figsize':(11,8)}
    #plotfigure.show = False
    xlimits = [-123.89,-123.73]
    ylimits = [46.94, 46.996]

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.axescmd = 'axes([.1,.1,.85,.5])'
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits


    def fixup(current_data):
        import pylab
        #addgauges(current_data)
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=range(54,58), format_string='ko', add_labels=True)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Aberdeen at t = %.2f hours' % t, fontsize=15)
        pylab.xticks(rotation=20,fontsize=10)
        pylab.yticks(fontsize=10)
        pylab.gca().set_aspect(1/pylab.cos(47*pylab.pi/180))
        pylab.ticklabel_format(useOffset=False)
        plot_src_boxes()

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin_water
    plotitem.pcolor_cmax = cmax_water
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'extend':'both'}
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    plotaxes = plotfigure.new_plotaxes('1d_plot')
    plotaxes.axescmd = 'axes([.1,.75,.7,.15])'
    plotaxes.afteraxes = plot_tide

    #-----------------------------------------
    # Figure for zoom
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Aberdeen speed', figno=11)
    plotfigure.kwargs = {'figsize':(11,8)}
    #plotfigure.show = False


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.axescmd = 'axes([.1,.1,.85,.5])'
    plotaxes.title = 'Speed'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits

    def fixup(current_data):
        import pylab
        #addgauges(current_data)
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=range(54,58), format_string='ko', add_labels=True)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Aberdeen speed at t = %.2f hours' % t, fontsize=15)
        pylab.xticks(rotation=20,fontsize=10)
        pylab.yticks(fontsize=10)
        pylab.gca().set_aspect(1/pylab.cos(47*pylab.pi/180))
        pylab.ticklabel_format(useOffset=False)
        plot_src_boxes()
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = speed
    plotitem.pcolor_cmap = colormaps.blue_yellow_red
    plotitem.pcolor_cmin = cmin_speed
    plotitem.pcolor_cmax = cmax_speed
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'extend':'max'}
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    plotaxes = plotfigure.new_plotaxes('1d_plot')
    plotaxes.axescmd = 'axes([.1,.75,.7,.15])'
    plotaxes.afteraxes = plot_tide

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'



    #-----------------------------------------
    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):
        from clawpack.visclaw import plot_timing_stats
        import os,sys
        try:
            timing_plotdir = plotdata.plotdir + '/_timing_figures'
            os.system('mkdir -p %s' % timing_plotdir)
            # adjust units for plots based on problem:
            units = {'comptime':'seconds', 'simtime':'hours', 
                     'cell':'millions'}
            plot_timing_stats.make_plots(outdir=plotdata.outdir, 
                                          make_pngs=True,
                                          plotdir=timing_plotdir, 
                                          units=units)
        except:
            print('*** Error making timing plots')

    otherfigure = plotdata.new_otherfigure(name='timing plots',
                    fname='_timing_figures/timing.html')
    otherfigure.makefig = make_timing_plots


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                # make multiple frame png's at once

    return plotdata

