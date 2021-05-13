# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:39:13 2021

@author: Loïc James McKeever
"""
import random as rd
import numpy as np

from bokeh.layouts import column, row
from bokeh.models import CustomJS, Slider, HoverTool
from bokeh.plotting import figure, output_file, show, ColumnDataSource

def MLE_LS_curve_fitting(hist, edges):
    hist_err = [(item - np.sqrt(item),item + np.sqrt(item)) for item in hist]
    hist_err_val = [np.sqrt(item) for item in hist]
    
    plot = figure(title="Number of Muon Decays in relation to decay time",width=1500,height=800)
    plot.quad(top=hist, bottom=hist, left=edges[:-1], right=edges[1:])
    plot.xaxis.axis_label = "Time in microseconds"
    plot.yaxis.axis_label = "Number of decays"
    
    #Create the MLE curvefit based on tau including the slider and the value of ln(L)
    x = [(edges[j-1]+edges[j])/2 for j in range(1,401)]
    y = [((3000*0.05)/2.5)*np.exp((-x[j-1])/2.5) for j in range(1,401)]
       
    lnL=sum([(hist[k-1] * np.log(y[k-1]*0.05) - (y[k-1]*0.05)) for k in range(1,len(hist)+1)])
    
    #Add the errorbars to the plot
    x_err = [(x_item,x_item) for x_item in x]
    plot.multi_line(x_err, hist_err)
    
    source = ColumnDataSource(data=dict(x=x,y=y,hist=hist))
    
    lnL_source = ColumnDataSource(data=dict(x=[2.5],y=[lnL],color=['#ff0000'],size=[10]))
    lnL_plot = figure(title="Ln(L) in reltation to tau", width=400, height=400)
    lnL_plot.scatter('x','y',color='color', size='size', source=lnL_source)
    lnL_plot.xaxis.axis_label = "Value of tau"
    lnL_plot.yaxis.axis_label = "Value of ln(L)"
    
    lnL_hovertool = HoverTool(tooltips=[("ln(L)","@y"),("tau","@x")])
    lnL_plot.tools.append(lnL_hovertool)
    
    plot.line('x','y', source=source, line_width=2, line_color='#ff0000', legend_label='MLE')
    
    tau_slider = Slider(start=0, end=5, value=2.5, step=.01, title="Tau")
    
    callback = CustomJS(args=dict(source=source, lnL_source=lnL_source, tau=tau_slider), code = """
        const data = source.data;
        const t = tau.value;
        const x = data['x'];
        const y = data['y'];
        const hist = data['hist'];
        
        const lnL_data = lnL_source.data;
        var clr = lnL_data['color'];
        var sz = lnL_data['size'];
        var lnL_x = lnL_data['x'];
        var lnL_y = lnL_data['y'];
        var lnL = 0;
        
        for (var i = 0; i < x.length; i++){
                y[i] = (3000*0.05)/t * Math.exp((-x[i])/t);
        }
        
        for (var i = 0; i < x.length; i++){
                lnL += hist[i] * Math.log(y[i]*0.05) - (y[i]*0.05);
        }
        
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
                
        source.change.emit();
        lnL_source.change.emit();
    """)
    
    tau_slider.js_on_change('value',callback)
    
    #Create the LS curve fit
    lnbins = [np.log(item) if item > 0 else 0 for item in hist]
    y_ls = [np.log(item) if item > 0 else 0 for item in y]
    
    B = []
    for i in range(0,400):
        l = (lnbins[i]-y_ls[i])**2
        h = (hist_err_val[i]/hist[i])**2
        if str(h) == 'nan':
            B.append(0)
        else:
            B.append(l/h)
    
    chi2 = sum(B)
    
    source_ls = ColumnDataSource(data=dict(x=x,y=y,y_ls=y_ls,lnbins=lnbins,hist=hist,hist_err_val=hist_err_val,))
    
    chi2_source = ColumnDataSource(data=dict(x=[2.5],y=[chi2],color=['#ff0000'],size=[10]))
    chi2_plot = figure(title="X^2 in reltation to tau", width=400, height=400)
    chi2_plot.scatter('x','y',color='color', size='size', source=chi2_source)
    chi2_plot.xaxis.axis_label = "Value of tau"
    chi2_plot.yaxis.axis_label = "Value of X^2"
    
    chi2_hovertool = HoverTool(tooltips=[("X^2","@y"),("tau","@x")])
    chi2_plot.tools.append(chi2_hovertool)
    
    plot.line('x','y', source=source_ls, line_width=2, line_color='#ffa500', legend_label='LS')
    
    tau_slider_ls = Slider(start=0, end=5, value=2.5, step=.001, title="Tau")
    
    callback_ls = CustomJS(args=dict(source=source_ls, chi2_source=chi2_source, tau=tau_slider_ls), code = """
        const data = source.data;
        const t = tau.value;
        const x = data['x'];
        const y = data['y'];
        const hist = data['hist'];
        const y_ls = data['y_ls'];
        const lnbins = data['lnbins'];
        const hist_err_val = data['hist_err_val']
        
        const chi2_data = chi2_source.data;
        var clr = chi2_data['color'];
        var sz = chi2_data['size'];
        var chi2_x = chi2_data['x'];
        var chi2_y = chi2_data['y'];
        var B = [];
        var chi2 = 0;
        
        for (var i = 0; i < x.length; i++){
                y[i] = (3000*0.05)/t * Math.exp((-x[i])/t);
        }
        
        for (var i = 0; i < x.length; i++){
                if (y[i] > 0){
                        y_ls[i] = Math.log(y[i]);
                } else {
                        y_ls[i] = 0;
                }
        }
                
        for (var i = 0; i < x.length; i++){
                var l = (lnbins[i]-y_ls[i])**2;
                var h = (hist_err_val[i]/hist[i])**2;
                B.push(l/h);     
        }
        
        for (var i = 0; i < x.length; i++){
                if (!isNaN(B[i])){
                        chi2 += B[i];
                }
        }
        
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
        source.change.emit();
        chi2_source.change.emit();
    """)
    
    tau_slider_ls.js_on_change('value',callback_ls)
    
    #Legend location and interaction option
    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"
    
    return plot, tau_slider, lnL_plot, tau_slider_ls, chi2_plot

#Experimental data plotting
file = open("LevangieMcKeever_3000.data", "r")
file = file.read()

lines = file.split("\n")
lines.pop()

decays = [int(line.split(" ")[0])/1000 for line in lines if int(line.split(" ")[0]) < 40000]

hist, edges = np.histogram(decays,bins=400, range=(0,max(decays)))
plot, tau_slider, lnL_plot, tau_slider_ls, chi2_plot = MLE_LS_curve_fitting(hist, edges)

#Simulated data plotting
tau = 2.2
sim_decays = [-tau*np.log(rd.random()) for i in range(0,3000)]

sim_hist, sim_edges=np.histogram(sim_decays, bins=400, range=(0,max(sim_decays)))
sim_plot, sim_tau_slider, sim_lnL_plot, sim_tau_slider_ls, sim_chi2_plot = MLE_LS_curve_fitting(sim_hist, sim_edges)

# sim_plot = figure(title="Simulated Muon Decays",width=1500,height=800)
# sim_plot.quad(top=sim_hist, bottom=0, left=sim_edges[:-1], right=sim_edges[1:],fill_color='#58e07c')

output_file("Muon_Decay.html", title="Muon Decay")

layout = column(row(plot, column(tau_slider, lnL_plot, tau_slider_ls, chi2_plot)), 
                row(sim_plot, column(sim_tau_slider, sim_lnL_plot, sim_tau_slider_ls, sim_chi2_plot)))
show(layout)