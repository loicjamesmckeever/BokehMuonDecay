# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:39:13 2021

@author: Loïc James McKeever
"""
import random as rd
import numpy as np
import sys
#Output errors to py.log file
sys.stderr = open('py.log', 'w')

from bokeh.layouts import column, row
from bokeh.models import CustomJS, Slider, HoverTool, Label, Button, Div, Range1d
from bokeh.plotting import figure, output_file, show, ColumnDataSource

def calculate_chi2(lnbins, y_ls_tau, hist_err_val, hist):
    B = []
    for i in range(0,400):
        l = (lnbins[i]-y_ls_tau[i])**2
        if hist[i] == 0:
            B.append(0)
        else:
            h = (hist_err_val[i]/hist[i])**2
            B.append(l/h)
    chi2 = sum(B)
    return chi2

def MLE_LS_curve_fitting(hist, edges):
    #Calculate the error for each value in the histogram
    hist_err = [(item - np.sqrt(item),item + np.sqrt(item)) for item in hist]
    hist_err_val = [np.sqrt(item) for item in hist]

    plot = figure(title="Number of Muon Decays in relation to decay time",width=1500,height=800)
    plot.quad(top=hist, bottom=hist, left=edges[:-1], right=edges[1:])
    plot.xaxis.axis_label = "Time in microseconds"
    plot.yaxis.axis_label = "Number of decays"

    #Create the MLE curvefit based on tau including the slider and the value of ln(L)
    #x is the time, each value is the center of the edges of the bins of the histogram
    #need a list of 500 of them for the lists being sent to ColumnDataSource to match
    x = [(edges[j-1]+edges[j])/2 for j in range(1,401)]
    x_datasource = [x for i in range(0,500)]
    #create the tau slider
    tau_slider = Slider(start=0.01, end=5, value=2.5, step=.01, title="Tau")
    #y is our MLE curve fit for each value of tau of our slider
    #Needs to be in list format for Bokeh
    slider_values = np.arange(0.01, 5.01, .01)
    y = [[((3000*0.05)/slider_value)*np.exp((-x[j-1])/slider_value) for j in range(1,401)] for slider_value in slider_values]
    #Calculate lnL for each value of tau between .01 and 5
    lnL_nans = [sum([(hist[k-1] * np.log(y[m][k-1]*0.05) - (y[m][k-1]*0.05)) for k in range(1,len(hist)+1)]) for m in range(0,len(slider_values))]
    #Get rid of any nans in lnL to keep JS happy
    lnL = [l if str(l) != 'nan' else 0 for l in lnL_nans]

    #Add the errorbars to the plot
    x_err = [(x_item,x_item) for x_item in x]
    plot.multi_line(x_err, hist_err)

    #Prepare the data for use with JS in the HTML file
    source = ColumnDataSource(data=dict(x=x_datasource,y=y, lnL=lnL))
    plot_source = ColumnDataSource(data=dict(x=x,y=y[250]))

    lnL_source = ColumnDataSource(data=dict(x=[slider_values[250]],y=[lnL[250]],color=['#ff0000'],size=[10]))
    lnL_plot = figure(title="Ln(L) in relation to tau", width=400, height=400)
    lnL_plot.scatter('x','y',color='color', size='size', source=lnL_source)
    lnL_plot.xaxis.axis_label = "Value of tau"
    #lnL_plot.x_range = Range1d(0, 5)
    lnL_plot.yaxis.axis_label = "Value of ln(L)"
    #lnL_plot.y_range = Range1d(0, 1000)

    lnL_hovertool = HoverTool(tooltips=[("ln(L)","@y"),("tau","@x")])
    lnL_plot.tools.append(lnL_hovertool)

    plot.line('x','y', source=plot_source, line_width=2, line_color='#ff0000', legend_label='MLE')

    lnL_label = Label(x=70, y=70, x_units='screen', y_units='screen',
                 text='ln(L) = ' + str(round(lnL[250], 2)), render_mode='css',
                 border_line_color='black', border_line_alpha=1.0,
                 background_fill_color='white', background_fill_alpha=1.0)

    #JavaScript to control plot based on slider action
    callback = CustomJS(args=dict(source=source, plot_source=plot_source, lnL_source=lnL_source, tau=tau_slider, lnL_label=lnL_label), code = """
        const data = source.data;
        const x_array = data['x'];
        const y_array = data['y'];
        const lnL_array = data['lnL']

        const t = tau.value;

        const plot_data = plot_source.data
        const x = plot_data['x']

        const lnL_data = lnL_source.data;
        var clr = lnL_data['color'];
        var sz = lnL_data['size'];
        var lnL_x = lnL_data['x'];
        var lnL_y = lnL_data['y'];
        var lnL = lnL_data['y'];

        var index = (t*100).toFixed()

        lnL = lnL_array[index]
        plot_data['y'] = y_array[index]

        lnL_label.text = 'ln(L) = ' + lnL.toFixed(2);

        if (!isNaN(lnL) && !lnL_x.includes(t) ) {
                lnL_y.push(lnL);
                lnL_x.push(t);
        }

        for (var i = 0; i < lnL_y.length; i++){
                if (lnL_y[i] == Math.max(...lnL_y)){
                        clr[i] = '#ff0000';
                        sz[i] = 10;
                } else {
                        clr[i] = '#0000ff';
                        sz[i] = 4;
                }
        }

        plot_source.change.emit();
        lnL_source.change.emit();
    """)

    tau_slider.js_on_change('value',callback)

    button = Button(label="Maximize ln(L)", button_type="success")
    button.js_on_click(CustomJS(args=dict(tau=tau_slider, source=source), code="""
        const data = source.data;
        const lnL_array = data['lnL']

        const max = Math.max(...lnL_array)
        const index = (lnL_array.indexOf(max)/100).toFixed(2)
        while (tau.value != index){
            if (tau.value < index){
                var tau_temp = tau.value;
                tau_temp = Number(tau_temp)
                tau_temp += .01
                tau.value = tau_temp.toFixed(2)
            } else if (tau.value > index){
                var tau_temp = tau.value;
                tau_temp = Number(tau_temp)
                tau_temp -= .01
                tau.value = tau_temp.toFixed(2)
            }
        }
    """))

    lnL_plot.add_layout(lnL_label)

    #Create the LS curve fit
    #List of natural log of each value in the histogram
    lnbins = [np.log(item) if item > 0 else 0 for item in hist]
    #List of natural log of each item of the curve fit for each value of tau
    y_ls = [[np.log(item) if item > 0 else 0 for item in y_tau] for y_tau in y]
    #Calculate chi2 for each value of tau
    chi2s = [calculate_chi2(lnbins, y_ls_tau, hist_err_val, hist) for y_ls_tau in y_ls]

    source_ls = ColumnDataSource(data=dict(x=x_datasource,y=y,chi2s=chi2s))
    plot_source_ls = ColumnDataSource(data=dict(x=x,y=y[250]))

    chi2_source = ColumnDataSource(data=dict(x=[2.5],y=[chi2s[250]],color=['#ff0000'],size=[10]))
    chi2_plot = figure(title="X^2 in relation to tau", width=400, height=400)
    chi2_plot.scatter('x','y',color='color', size='size', source=chi2_source)
    chi2_plot.xaxis.axis_label = "Value of tau"
    chi2_plot.yaxis.axis_label = "Value of X^2"

    chi2_hovertool = HoverTool(tooltips=[("X^2","@y"),("tau","@x")])
    chi2_plot.tools.append(chi2_hovertool)

    plot.line('x','y', source=plot_source_ls, line_width=2, line_color='#ffa500', legend_label='LS')

    tau_slider_ls = Slider(start=0.01, end=5, value=2.5, step=.01, title="Tau")

    chi2_label = Label(x=70, y=70, x_units='screen', y_units='screen',
                 text='X^2 = ' + str(round(chi2s[250], 2)), render_mode='css',
                 border_line_color='black', border_line_alpha=1.0,
                 background_fill_color='white', background_fill_alpha=1.0)

    callback_ls = CustomJS(args=dict(source=source_ls, plot_source=plot_source_ls, chi2_source=chi2_source, tau=tau_slider_ls, chi2_label=chi2_label), code = """
        const data = source.data;
        const x_array = data['x'];
        const y_array = data['y'];
        const chi2_array = data['chi2s'];

        const t = tau.value;

        const plot_data = plot_source.data
        const x = plot_data['x']

        const chi2_data = chi2_source.data;
        var clr = chi2_data['color'];
        var sz = chi2_data['size'];
        var chi2_x = chi2_data['x'];
        var chi2_y = chi2_data['y'];
        var chi2 = chi2_data['y'];

        var index = (t*100).toFixed()

        chi2 = chi2_array[index]
        plot_data['y'] = y_array[index]

        chi2_label.text = 'X^2 = ' + chi2.toFixed(2);

        if (!isNaN(chi2) && !chi2_x.includes(t) ) {
                chi2_y.push(chi2);
                chi2_x.push(t);
        }

        for (var i = 0; i < chi2_y.length; i++){
                if (chi2_y[i] == Math.min(...chi2_y)){
                        clr[i] = '#ff0000';
                        sz[i] = 10;
                } else {
                        clr[i] = '#0000ff';
                        sz[i] = 4;
                }
        }
        plot_source.change.emit();
        chi2_source.change.emit();
    """)

    tau_slider_ls.js_on_change('value',callback_ls)

    chi2_button = Button(label="Minimize X^2", button_type="success")
    chi2_button.js_on_click(CustomJS(args=dict(tau=tau_slider_ls, source=source_ls), code="""
        const data = source.data;
        const chi2_array = data['chi2s']

        const min = Math.min(...chi2_array)
        const index = (chi2_array.indexOf(min)/100).toFixed(2)

        while (tau.value != index){
            if (tau.value < index){
                var tau_temp = tau.value;
                tau_temp = Number(tau_temp)
                tau_temp += .01
                tau.value = tau_temp.toFixed(2)
            } else if (tau.value > index){
                var tau_temp = tau.value;
                tau_temp = Number(tau_temp)
                tau_temp -= .01
                tau.value = tau_temp.toFixed(2)
            }
        }
    """))

    chi2_plot.add_layout(chi2_label)

    #Legend location and interaction option
    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"

    return plot, tau_slider, lnL_plot, button, tau_slider_ls, chi2_plot, chi2_button

#Experimental data plotting
file = open("LevangieMcKeever_3000.data", "r")
file = file.read()

lines = file.split("\n")
lines.pop()

decays = [int(line.split(" ")[0])/1000 for line in lines if int(line.split(" ")[0]) < 40000]

hist, edges = np.histogram(decays,bins=400, range=(0,20))
plot, tau_slider, lnL_plot, button, tau_slider_ls, chi2_plot, chi2_button = MLE_LS_curve_fitting(hist, edges)

#Simulated data plotting
tau = 2.2
sim_decays = [-tau*np.log(rd.random()) for i in range(0,3000)]

sim_hist, sim_edges=np.histogram(sim_decays, bins=400, range=(0,20))
sim_plot, sim_tau_slider, sim_lnL_plot, sim_button, sim_tau_slider_ls, sim_chi2_plot, sim_chi2_button = MLE_LS_curve_fitting(sim_hist, sim_edges)

#Define CSS style for HTML divs
style_title = {'font-size': '300%', 'font-family':'Georgia, serif'}
style_div = {'font-size': '200%', 'font-family':'Georgia, serif', 'width':'1900px'}
style_equation = {'font-size': '200%', 'text-align':'center', 'width':'1500px'}

#Various divs and headers for HTML output_file
header = Div(text = """
             <header> Muon Decay </header>
             """
             , style = style_title)

intro_header = Div(text = """
             <header> Introduction </header>
             """
             , style = style_div)

intro = Div(text = """
         <div>
         In this experiment we attempt to determine the mean lifetime of muons. This
         will be done by detecting muons and muon decays using a scintillator detector
         and specialized software. The data will be filtered so as to only use data points
         from the detection of decays. Once filtered the data is binned and curve fit using the Maximum
         Likelyhood and Least Squares methods. We will test both methods on a set of
         data we simulate ourselves using the known mean lifetime as well as the data set from our detector.
         </div>
         """
         , style = style_div)

#Setting up output file and layout
output_file("Muon_Decay_PY.html", title="Muon Decay")

layout = column(header, intro_header, intro, row(plot, column(tau_slider, lnL_plot, button, tau_slider_ls, chi2_plot, chi2_button)),
                row(sim_plot, column(sim_tau_slider, sim_lnL_plot, sim_button, sim_tau_slider_ls, sim_chi2_plot, sim_chi2_button)))

show(layout)
