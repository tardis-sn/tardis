from bqplot import *
from bqplot.interacts import (BrushSelector)
from bqplot import LinearScale, Axis, Lines, Scatter, Figure
from ipywidgets import VBox,HTML
from bqplot import pyplot as plt
from ipywidgets import widgets
from ipywidgets import *
from ipywidgets import Layout

def line_Plot(mdl):
    mdl=mdl
    x_sc = LinearScale()
    y_sc = LinearScale()


    x_data = mdl.runner.spectrum.wavelength
    y_data = mdl.runner.spectrum.luminosity_density_lambda
    z_data = mdl.runner.spectrum_virtual.luminosity_density_lambda
    a_data = mdl.runner.spectrum_integrated.luminosity_density_lambda
    
    
    symbol = 'Spectrum Wavelenght'
    symbol2 = 'Spectrum Density'

    #def_tt = Tooltip(fields=['x_data', 'y_data'], formats=['', '.4f'])
    line_chart1 = Lines(x=x_data, y=y_data, scales= {'x': x_sc, 'y': y_sc}, 
                            display_legend=True, labels=["Normal Packets"] )

    line_chart2 = Lines(x=x_data, y=z_data, scales= {'x': x_sc, 'y': y_sc}, 
                            display_legend=True, labels=["Virtual Packets"] , colors=['red'])
    line_chart3 = Lines(x=x_data, y=a_data, scales= {'x': x_sc, 'y': y_sc}, 
                            display_legend=True, labels=["Formal Integral"] , colors=['green'])

    ax_x = Axis(scale=x_sc,label="Wavelength [$\AA$]")
    ax_y = Axis(scale=y_sc, label='Luminosity [erg/s/$\AA$]',orientation='vertical')

    br_sel = BrushSelector(x_scale=x_sc, y_scale=y_sc, marks=[line_chart1 ,line_chart2, line_chart3], color='white')
    
    db_scat_brush = HTML(value='[]')
    
    
    ## Now, we define a function that will be called when the selectors are interacted with - a callback
    ## call back for the selector
    def brush_callback(change):
        db_scat_brush.value = 'The selected Values are ' + str(br_sel.selected)
    
    br_sel.observe(brush_callback, names=['brushing'])
    
    panzoom = PanZoom(scales= {'x': [x_sc] , 'y': [y_sc]})
    
    fig_scat_brush3 =Figure(marks=[line_chart1 , line_chart2 , line_chart3], axes=[ax_x, ax_y],title='Brush Selector',
                            interaction=br_sel,background_style={'fill':'moccasin'})
    
    
    from collections import OrderedDict
    from traitlets import link
    
    toggle = widgets.ToggleButtons(options=OrderedDict([('BrushIntervalSelector', br_sel),('Zoom', panzoom),('None', None)]))
    
    link((toggle, 'value'), (fig_scat_brush3, 'interaction'))
    return VBox([db_scat_brush,fig_scat_brush3 ,toggle], align_self='stretch')    

def line_Plot_show(mdl):
    x_sc = LinearScale()
    y_sc = LinearScale()


    x_data = mdl.runner.spectrum.wavelength
    y_data = mdl.runner.spectrum.luminosity_density_lambda
    z_data = mdl.runner.spectrum_virtual.luminosity_density_lambda
    a_data = mdl.runner.spectrum_integrated.luminosity_density_lambda

    #def_tt = Tooltip(fields=['x_data', 'y_data'], formats=['', '.4f'])
    line_chart1 = Lines(x=x_data, y=y_data, scales= {'x': x_sc, 'y': y_sc}, 
                            display_legend=True, labels=["Normal Packets"] )

    line_chart2 = Lines(x=x_data, y=z_data, scales= {'x': x_sc, 'y': y_sc}, 
                            display_legend=True, labels=["Virtual Packets"] , colors=['red'])
    line_chart3 = Lines(x=x_data, y=a_data, scales= {'x': x_sc, 'y': y_sc}, 
                            display_legend=True, labels=["Formal Integral"] , colors=['green'])

    ax_x = Axis(scale=x_sc)
    ax_y = Axis(scale=y_sc, orientation='vertical')

    plt.figure(marks=[line_chart1 , line_chart2 , line_chart3], axes=[ax_x, ax_y])
    plt.show()
    
def m_details(mdl):
    print("dilution factor in each cell : ",mdl.iterations_w)  # dilution factor in each cell
    print("radiation temperature in each cell : ",mdl.iterations_t_rad)  # radiation temperature in each cell
    print("electron density in each cell : ",mdl.iterations_electron_densities)  # electron density in each cell
    print("inner boundary temperature : ",mdl.iterations_t_inner)  # inner boundary temperature
    
def model_parameters(mdl):
    print("Iterations Executed : ",mdl.iterations)
    print("Model Converged : ",mdl.converged)
    print("Number of Packets : ",mdl.no_of_packets)
    print("Number of Virtual Packets : ",mdl.no_of_virtual_packets)
    print("Numbers of Last Number of Packets : ",mdl.last_no_of_packets)
    print("Inner Luminosity : ",mdl.luminosity_requested)
    print("Inner Temperature : ",mdl.iterations_t_inner.value.max()," K")
def plot_w(model):
    
    x_scc = LinearScale()
    y_scc = LinearScale()
    
    a=range(len(model.model.t_rad.value))
    b=model.model.t_rad
    x_data1 = a
    y_data2 = b


    line_chart2 = Lines(x=x_data1, y=y_data2, scales= {'x': x_scc, 'y': y_scc}, 
                            display_legend=True, labels=["T_rad"] ,colors=['red'])

    ax_x = Axis(label="Shell",scale=x_scc)
    ax_y = Axis(label="Rad_Temp",scale=y_scc, orientation='vertical')

    plt.figure(marks=[line_chart2], axes=[ax_x, ax_y],title='RAD Temp Vs Shell')
    return plt.show()

import pandas as pd
import numpy as np
from IPython.display import display

def table_data(mdl):
    data=mdl.iterations_w
    data = pd.DataFrame(data=data[0:,0:],
               index=[i for i in range(data.shape[0])],
               columns=['Iter_'+str(i) for i in range(data.shape[1])])

    return (display('Dilution Factor in each cell',data))

def table_data1(mdl):
    #print("radiation temperature in each cell : ",sim.iterations_t_rad)  # radiation temperature in each cell
    #print(type(sim.iterations_t_rad))
    data1=mdl.iterations_t_rad
    data1 = pd.DataFrame(data=data1[0:,0:],
                   index=[i for i in range(data1.shape[0])],
                   columns=['Rad_Temp_'+str(i) for i in range(data1.shape[1])])

    return (display('Radiation Temperature In Each Cell',data1))

def table_data2(mdl):
    #print(type(sim.iterations_electron_densities))
    data2=mdl.iterations_electron_densities
    data2 = pd.DataFrame(data=data2[0:,0:],
                   index=[i for i in range(data2.shape[0])],
                   columns=['Elec_Dens_'+str(i) for i in range(data2.shape[1])])

    return (display('Electron Density In Each Cell',data2))

def table_data3(mdl):
    #print("inner boundary temperature : ",sim.iterations_t_inner)  # inner boundary temperature
    #print(type(sim.iterations_t_inner))
    a=np.array(mdl.iterations_t_inner)
    data3=a.ravel()
    #print(data3)
    #print(type(data3))
    #print(data3.shape)
    #print(data3[:,])
    #data3=pd.Series(data3)
    data3=pd.Series(data3,name='Inner boundary temp')
    #print(data3)
    data3=pd.DataFrame(data3)
    return (display('Inner Boundary Temperature',data3))

def table_data(mdl):
    data=mdl.iterations_w
    data = pd.DataFrame(data=data[0:,0:],
               index=[i for i in range(data.shape[0])],
               columns=['Iter_'+str(i) for i in range(data.shape[1])])

    return (display('Dilution Factor in each cell',data))

def table_data_all(mdl):

    data=mdl.iterations_w
    data = pd.DataFrame(data=data[0:,0:],
               index=[i for i in range(data.shape[0])],
               columns=['Iter_'+str(i) for i in range(data.shape[1])])

    #print("radiation temperature in each cell : ",sim.iterations_t_rad)  # radiation temperature in each cell
    #print(type(sim.iterations_t_rad))
    data1=mdl.iterations_t_rad
    data1 = pd.DataFrame(data=data1[0:,0:],
                   index=[i for i in range(data1.shape[0])],
                   columns=['Rad_Temp_'+str(i) for i in range(data1.shape[1])])


    #print(type(sim.iterations_electron_densities))
    data2=mdl.iterations_electron_densities
    data2 = pd.DataFrame(data=data2[0:,0:],
                   index=[i for i in range(data2.shape[0])],
                   columns=['Elec_Dens_'+str(i) for i in range(data2.shape[1])])

    #print("inner boundary temperature : ",sim.iterations_t_inner)  # inner boundary temperature
    #print(type(sim.iterations_t_inner))
    a=np.array(mdl.iterations_t_inner)
    data3=a.ravel()
    #print(data3)
    #print(type(data3))
    #print(data3.shape)
    #print(data3[:,])
    #data3=pd.Series(data3)
    data3=pd.Series(data3,name='Inner boundary temp')
    #print(data3)
    data3=pd.DataFrame(data3)

    all_in_one = pd.concat([data,data1,data2,data3],axis=1)
    #all_in_one.head()

    return display("Access information about individual iterations",all_in_one)

def heat_plot_rad(mdl):
    data1=mdl.iterations_t_rad
    data1 = pd.DataFrame(data=data1[0:,0:],
                   index=[i for i in range(data1.shape[0])],
                   columns=['Rad_Temp_'+str(i) for i in range(data1.shape[1])])
    x = data1
    plt.figure(title='Heat_Map of Radiation temperature in each cell',padding_y=0)
    axes_options = {'x': {'label': 'T(K)'}, 'y': {'label':'Shell No'}}
    plt.heatmap(x,axes_options=axes_options)
    return plt.show()

def heat_plot_elec(mdl):

    data2=mdl.iterations_electron_densities
    data2 = pd.DataFrame(data=data2[0:,0:],
                   index=[i for i in range(data2.shape[0])],
                   columns=['Elec_Dens_'+str(i) for i in range(data2.shape[1])])
    y = data2
    plt.figure(title="Heat Map of Electron Density in each cell",padding_y=0)
    axes_options = {'x': {'label': 'T(K)'}, 'y': {'label':'Shell No'}}
    plt.heatmap(y,axes_options=axes_options)
    return plt.show()
