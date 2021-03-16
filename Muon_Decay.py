# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:39:13 2021

@author: Loïc James McKeever
"""
import random as rd
import numpy as np

from bokeh.layouts import column, row
from bokeh.models import CustomJS, Slider
from bokeh.plotting import figure, output_file, show, ColumnDataSource

#Experimental data plotting
file = open("Downloads\LevangieMcKeever_3000.data", "r")
file = file.read()

lines = file.split("\n")
lines.pop()

decays = [int(line.split(" ")[0])/1000 for line in lines if int(line.split(" ")[0]) < 40000]

hist, edges = np.histogram(decays,bins=400)

plot = figure(title="Muon Decays",width=1500,height=800)
plot.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])

#Create the curvefit based on tau including the slider and the value of ln(L)
xL = np.linspace(0,5,10000)
n = 5/10000

x = []
y = []
for j in range(1,401):
    x.append((0.025 + j*0.05))
    y.append(((3000*0.05)/2.5)*np.exp((-x[j-1])/2.5))
    
lnL=sum([(hist[k-1] * np.log(y[k-1]*0.05) - (y[k-1]*0.05)) for k in range(1,len(hist)+1)])

source = ColumnDataSource(data=dict(x=x,y=y,hist=hist))

lnL_source = ColumnDataSource(data=dict(x=[],y=[]))
lnL_plot = figure(title="Ln(L) in reltation to tau", width=400, height=400)
lnL_plot.scatter('x','y',source=lnL_source)

plot.line('x','y', source=source, line_width=2, line_color='#ff0000')

tau_slider = Slider(start=0, end=5, value=2.5, step=.01, title="Tau")

callback = CustomJS(args=dict(source=source, lnL_source=lnL_source, tau=tau_slider), code = """
    const data = source.data;
    const t = tau.value;
    const x = data['x'];
    const y = data['y'];
    const hist = data['hist'];
    
    const lnL_data = lnL_source.data;
    var lnL_x = lnL_data['x'];
    var lnL_y = lnL_data['y'];
    var lnL = 0;
    
    for (var i = 0; i < x.length; i++){
            y[i] = (3000*0.05)/t * Math.exp((-x[i])/t)
    }
    
    for (var i = 0; i < x.length; i++){
            lnL += hist[i] * Math.log(y[i]*0.05) - (y[i]*0.05);
    }
    
    lnL_y.push(lnL);
    lnL_x.push(t);
    
    console.log(lnL_x,lnL_y);
    
    source.change.emit();
    lnL_source.change.emit();
""")

tau_slider.js_on_change('value',callback)

#Simulated data plotting
tau = 2200
sim_decays = [-tau*np.log(rd.random()) for i in range(0,3000)]

sim_hist, sim_edges=np.histogram(sim_decays, bins=1000)

sim_plot = figure(title="Simulated Muon Decays",width=1500,height=800)
sim_plot.quad(top=sim_hist, bottom=0, left=sim_edges[:-1], right=sim_edges[1:],fill_color='#58e07c')

output_file("Muon_Decay.html", title="Muon Decay")

layout = column(row(plot, column(tau_slider, lnL_plot)), sim_plot)
show(layout)