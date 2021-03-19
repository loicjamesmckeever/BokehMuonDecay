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

#Experimental data plotting
file = open("LevangieMcKeever_3000.data", "r")
file = file.read()

lines = file.split("\n")
lines.pop()

decays = [int(line.split(" ")[0])/1000 for line in lines if int(line.split(" ")[0]) < 40000]

hist, edges = np.histogram(decays,bins=400, range=(0,max(decays)))

hist_err = [(item - np.sqrt(item),item + np.sqrt(item)) for item in hist]

plot = figure(title="Muon Decays",width=1500,height=800)
plot.quad(top=hist, bottom=hist, left=edges[:-1], right=edges[1:])

#Create the MLE curvefit based on tau including the slider and the value of ln(L)
x = []
y = []
for j in range(1,401):
    x.append((edges[j-1]+edges[j])/2)
    y.append(((3000*0.05)/2.5)*np.exp((-x[j-1])/2.5))
   
lnL=sum([(hist[k-1] * np.log(y[k-1]*0.05) - (y[k-1]*0.05)) for k in range(1,len(hist)+1)])

#Add the errorbars to the plot
x_err = [(x_item,x_item) for x_item in x]
plot.multi_line(x_err, hist_err)

source = ColumnDataSource(data=dict(x=x,y=y,hist=hist))

lnL_source = ColumnDataSource(data=dict(x=[],y=[],color=[],size=[]))
lnL_plot = figure(title="Ln(L) in reltation to tau", width=400, height=400)
lnL_plot.scatter('x','y',color='color', size='size', source=lnL_source)

lnL_hovertool = HoverTool(tooltips=[("ln(L)","@y"),("tau","@x")])
lnL_plot.tools.append(lnL_hovertool)

plot.line('x','y', source=source, line_width=2, line_color='#ff0000')

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



#Simulated data plotting
tau = 2200
sim_decays = [-tau*np.log(rd.random()) for i in range(0,3000)]

sim_hist, sim_edges=np.histogram(sim_decays, bins=1000)

sim_plot = figure(title="Simulated Muon Decays",width=1500,height=800)
sim_plot.quad(top=sim_hist, bottom=0, left=sim_edges[:-1], right=sim_edges[1:],fill_color='#58e07c')

output_file("Muon_Decay.html", title="Muon Decay")

layout = column(row(plot, column(tau_slider, lnL_plot)), sim_plot)
show(layout)