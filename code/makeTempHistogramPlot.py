from numpy import *
import math
import scipy 
import numpy as np
from pylab import *
import matplotlib.pyplot as pl
from matplotlib import rc
from matplotlib.ticker import NullFormatter

def ellipse(ra,rb,ang,xcenter,ycenter):
    theta_bins = linspace(start = ang, stop = 2.0*pi, num = 100)
    X = ra*cos(theta_bins) + xcenter
    Y = rb*sin(theta_bins) + ycenter
    return X,Y

def MakeTempHistogramPlot(xdata,ydata,filename=None,clues=[], clues_label='',xlims=-99, ylims =-99 , \
                              nxbins = 50,nybins=50, bw=0, nbins=100,contours=1,sigma=1,line=1, 
                          xlabel='A', ylabel='B', mean_obs=0, sigma_obs=0, 
                          X_field=array([]), Y_field=array([]), Z_field=array([]), levels=array([])):
#bw = 0 for color, = 1 for black and white
#line = 0 for no line, =1 for line
#sigma = 1 for display % below line, =0 for not
#contours = 1 for display 1,2,3 sigma contours, = 0 for not.
    
    

    # Define the x and y data
    x = xdata
    y = ydata

# Set up default x and y limits
    if (xlims == -99): xlims = [min(x),max(x)]
    if (ylims == -99): ylims = [min(y),max(y)]
    
# Set up your x and y labels
    xlabel = xlabel
    ylabel = ylabel
    mtitle = ''
    
    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    
    # Set up the geometry of the three plots
    rect_temperature = [left, bottom, width, height] # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
    
    # Set up the size of the figure
    fig = plt.figure(1, figsize=(9.5,9))
    
    # Make the three plots
    axTemperature = plt.axes(rect_temperature) # temperature plot
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram
    
    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Find the min/max of the data
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(ylims)

    # Make the 'main' temperature plot
    xbins = linspace(start = xmin, stop = xmax, num = nxbins)
    ybins = linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    aspectratio = 1.0*(xmax - xmin)/(1.0*ymax - ymin)
    H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins))
    X = xcenter
    Y = ycenter
    Z = H

    # Plot the temperature data
    if(bw): cax = axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax], \
                                           interpolation='nearest', origin='lower',aspect=aspectratio, cmap=cm.gist_yarg)
    else : cax = axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax], \
                                          interpolation='nearest', origin='lower',aspect=aspectratio)

    # Plot the temperature plot contours
    if(bw): contourcolor = 'black'
    else: contourcolor = 'white'

    if clues.size:
        axTemperature.scatter(clues[:,0], clues[:,1], color='white', edgecolors='black', s=200, label=clues_label)        
        axTemperature.legend(loc=0, scatterpoints=1, fontsize=20)

    if X_field.size:
        print 'HOLA'
        C=axTemperature.contour(X_field,Y_field,Z_field, levels, colors='k', linewidths=2)


    if(mean_obs!=0):
        X,Y=ellipse(sigma_obs[0], sigma_obs[1], 0, mean_obs[0], mean_obs[1])
        axTemperature.plot(X,Y,color = contourcolor,ms=1,linewidth=3.0)
        axTemperature.annotate('$\mathrm{Obs.}$', xy=(X[45], Y[45]), xycoords='data',xytext=(-10, 10), textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize=25, color = contourcolor)
    if (contours==0):
        print ''
    elif (contours==1):
        xcenter = np.mean(x)
        ycenter = np.mean(y)
        ra = np.std(x)
        rb = np.std(y)
        print xcenter, ycenter, ra, rb
        ang = 0
        X,Y=ellipse(ra,rb,ang,xcenter,ycenter)
        axTemperature.plot(X,Y,"k:",ms=1,linewidth=2.0)
        axTemperature.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10), textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize=25)
        X,Y=ellipse(2*ra,2*rb,ang,xcenter,ycenter)
        axTemperature.plot(X,Y,"k:",color = contourcolor,ms=1,linewidth=2.0)
        axTemperature.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10), textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize=25, color = contourcolor)
        X,Y=ellipse(3*ra,3*rb,ang,xcenter,ycenter)
        axTemperature.plot(X,Y,"k:",color = contourcolor, ms=1,linewidth=2.0)
        axTemperature.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10), textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize=25, color = contourcolor)
    else:
        xcenter = np.mean(x)
        ycenter = np.mean(y)
        ra = np.std(x)
        rb = np.std(y)
        ang = contours*np.pi/180.0
        X,Y=ellipse(ra,rb,ang,xcenter,ycenter)
        axTemperature.plot(X,Y,"k:",ms=1,linewidth=2.0)
        axTemperature.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10), textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize=25)
        X,Y=ellipse(2*ra,2*rb,ang,xcenter,ycenter)
        axTemperature.plot(X,Y,"k:",ms=1,linewidth=2.0, color = contourcolor)
        axTemperature.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10), textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize=25, color = contourcolor)
        X,Y=ellipse(3*ra,3*rb,ang,xcenter,ycenter)
        axTemperature.plot(X,Y,"k:",ms=1,linewidth=2.0, color = contourcolor)
        axTemperature.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data', xytext=(10, 10), textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize=25, color = contourcolor)

#Plot the % below line
#    belowline = 1.0*size(where((x - y) > 0.0))/size(x)*1.0*100
#    if(sigma): axTemperature.annotate('$%.2f\%%\mathrm{\\ Below\\ Line}$'%(belowline), xy=(xmax-100, ymin+3),fontsize=20, color = contourcolor)

#Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=25)
    axTemperature.set_ylabel(ylabel,fontsize=25)

#Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')

    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')

#Plot the line on the temperature plot
#    if(line): axTemperature.plot([-1000,1000], [-1000,1000], 'k-', linewidth=2.0, color = contourcolor)

#Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)

#Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)

#Plot the histograms
    if (bw):
        axHistx.hist(x, bins=xbins, color = 'silver')
        axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'dimgray')
    else:
        axHistx.hist(x, bins=xbins, color = 'blue')
        axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'red')

#Set up the histogram limits
    axHistx.set_xlim( min(xlims), max(xlims) )
    axHisty.set_ylim( min(ylims), max(ylims))

#Make the tickmarks pretty
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

#Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

#Cool trick that changes the number of tickmarks for the histogram axes
        axHisty.xaxis.set_major_locator(MaxNLocator(4))
        axHistx.yaxis.set_major_locator(MaxNLocator(4))

    if(filename):
        savefig(filename + '.eps',format = 'eps', transparent=True)
        savefig(filename + '.pdf',format = 'pdf', transparent=True)
        savefig(filename + '.png',format = 'png', transparent=True)
        
    return 0

