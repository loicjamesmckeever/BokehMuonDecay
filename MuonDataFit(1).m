%Get the file name from the user
filename=inputdlg('Data file to be analyzed?', 'File Selection');
filename=string(filename);
%Get 3000 muon decay data from file
fileID = fopen(filename);
[decaydata, nm] = fscanf(fileID,'%i');
fclose(fileID);
%Divide the data up into 50 ns bins
decaydata=decaydata./1000;
maxdecay=max(decaydata);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(decaydata,edges);
sz=size(edges);
count=sz(1,2);
bins=zeros(count-1,1);
i=1;
while i<count
    bins(i)=sum(databinned==i);
    i=i+1;
end
binerror=bins.^(1/2);
m=size(bins);
m=m(1,1);
j=1;
ploterror=zeros(m,1);
while j<m
    ploterror(j,1)=binerror(j,1);
    j=j+5;
end
%We assume the decay is somewhere between 0 and 5 microseconds
xL=linspace(0,5,10000)';
n=5/10000;
%Perform the ML to best fit the data with an exponential decay curve
%Determine ln(L) for tau between 0 and 5
lnL=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    expected=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        y(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        i=i+1;
    end
    lnL(j,1)=sum(bins.*log(y*0.05)-(y*0.05));
    j=j+1;
end
%Determine the highest value of ln(L) and set tau to that value
[maxlnL,maxindex]=max(lnL);
[~,max1index]=min(abs(lnL-(maxlnL-2)));
tau=xL(maxindex,1);
dtau=xL(max1index,1);
error=abs(tau-dtau);
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
figure (1)
set(figure(1),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.');
hold on
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('ML method curve fit of %i Decay muon decays in 50 ns bins.',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(1),'ML_Curve_Fit_3000.jpg');
%Least squares method fitting
lnbins=zeros(count-1,1);
i=1;
while i<count
    if sum(databinned==i)~=0
        lnbins(i)=log(sum(databinned==i));
    else
        lnbins(i)=0;
    end
    i=i+1;
end                        
m=size(lnbins);
m=m(1,1);
xL=linspace(0,5,10000)';
n=5/10000;
chi2=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    ys=zeros(m,1);
    B=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        ys(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        if (log(nm*0.05)-log(j*n)-((0.025+i*.05)/(j*n)))>0
            y(i,1)=log(nm*0.05)-log(j*n)-((0.025+i*.05)/(j*n));
        else
            y(i,1)=0;
        end
        i=i+1;
    end
    i=1;
    while i<m+1
        nan=isnan(((lnbins(i,1)-y(i,1))^2)/((binerror(i,1)/bins(i,1))^2));
        if nan==0
            B(i,1)=((lnbins(i,1)-y(i,1))^2)/((binerror(i,1)/bins(i,1))^2); 
        else
            B(i,1)=0;
        end
        i=i+1;
    end
    chi2(j,1)=sum(B);
    j=j+1;
end
%Determine the lowest value of chi2 and set tau to that value
[minchi2,minindex]=min(chi2);
[min1,min1index]=min(abs(chi2-(minchi2+1)));
tau=xL(minindex,1);
dtau=xL(min1index,1);
error=(abs(tau-dtau));
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
%edges(end)=[];
figure (2)
set(figure(2),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.');
hold on
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('LS method curve fit of %i Decay muon decays in 50 ns bins.',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(2),'LS_Curve_Fit_3000.jpg');
%Divide the data up into chunks of 500 points
nm=500;
i=1;
decaydata1=zeros(500,1);
decaydata2=zeros(500,1);
decaydata3=zeros(500,1);
decaydata4=zeros(500,1);
decaydata5=zeros(500,1);
decaydata6=zeros(500,1);
while i<500
    decaydata1(i,1)=decaydata(i,1);
    decaydata2(i,1)=decaydata(i+500,1);
    decaydata3(i,1)=decaydata(i+1000,1);
    decaydata4(i,1)=decaydata(i+1500,1);
    decaydata5(i,1)=decaydata(i+2000,1);
    decaydata6(i,1)=decaydata(i+2500,1);
    i=i+1;
end
%Perform a ML fit on each data set
%First Set
%Divide the data up into 50 ns bins
maxdecay=max(decaydata1);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(decaydata1,edges);
sz=size(edges);
count=sz(1,2);
bins=zeros(count-1,1);
i=1;
while i<count
    bins(i)=sum(databinned==i);
    i=i+1;
end
binerror=bins.^(1/2);
m=size(bins);
m=m(1,1);
ploterror=zeros(m,1);
j=1;
while j<m
    ploterror(j,1)=binerror(j,1);
    j=j+5;
end
%We assume the decay is somewhere between 0 and 5 microseconds
xL=linspace(0,5,10000)';
n=5/10000;
%Perform the ML to best fit the data with an exponential decay curve
%Determine ln(L) for tau between 0 and 5
lnL=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    expected=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        y(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        i=i+1;
    end
    lnL(j,1)=sum(bins.*log(y*0.05)-(y*0.05));
    j=j+1;
end
%Determine the highest value of ln(L) and set tau to that value
[~,max1index]=min(abs(lnL-(maxlnL-2)));
tau=xL(maxindex,1);
dtau=xL(max1index,1);
error=abs(tau-dtau);
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
figure (3)
set(figure(3),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.'); 
hold on 
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau, error);
tcaption=sprintf('ML method curve fit of %i Decay muon decays in 50 ns bins, first set of 500',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(3),'ML_Curve_Fit_500a.jpg');
%Second Set
%Divide the data up into 50 ns bins
maxdecay=max(decaydata2);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(decaydata2,edges);
sz=size(edges);
count=sz(1,2);
bins=zeros(count-1,1);
i=1;
while i<count
    bins(i)=sum(databinned==i);
    i=i+1;
end
binerror=bins.^(1/2);
m=size(bins);
m=m(1,1);
ploterror=zeros(m,1);
j=1;
while j<m
    ploterror(j,1)=binerror(j,1);
    j=j+5;
end
%We assume the decay is somewhere between 0 and 5 microseconds
xL=linspace(0,5,10000)';
n=5/10000;
%Perform the ML to best fit the data with an exponential decay curve
%Determine ln(L) for tau between 0 and 5
lnL=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    expected=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        y(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        i=i+1;
    end
    lnL(j,1)=sum(bins.*log(y*0.05)-(y*0.05));
    j=j+1;
end
%Determine the highest value of ln(L) and set tau to that value
[~,max1index]=min(abs(lnL-(maxlnL-2)));
tau=xL(maxindex,1);
dtau=xL(max1index,1);
error=abs(tau-dtau);
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
figure (4)
set(figure(4),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.'); 
hold on 
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('ML method curve fit of %i Decay muon decays in 50 ns bins, second set of 500',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(4),'ML_Curve_Fit_500b.jpg');
%Third Set
%Divide the data up into 50 ns bins
maxdecay=max(decaydata3);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(decaydata3,edges);
sz=size(edges);
count=sz(1,2);
bins=zeros(count-1,1);
i=1;
while i<count
    bins(i)=sum(databinned==i);
    i=i+1;
end
binerror=bins.^(1/2);
m=size(bins);
m=m(1,1);
ploterror=zeros(m,1);
j=1;
while j<m
    ploterror(j,1)=binerror(j,1);
    j=j+5;
end
%We assume the decay is somewhere between 0 and 5 microseconds
xL=linspace(0,5,10000)';
n=5/10000;
%Perform the ML to best fit the data with an exponential decay curve
%Determine ln(L) for tau between 0 and 5
lnL=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    expected=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        y(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        i=i+1;
    end
    lnL(j,1)=sum(bins.*log(y*0.05)-(y*0.05));
    j=j+1;
end
%Determine the highest value of ln(L) and set tau to that value
[~,max1index]=min(abs(lnL-(maxlnL-2)));
tau=xL(maxindex,1);
dtau=xL(max1index,1);
error=abs(tau-dtau);
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
figure (5)
set(figure(5),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.'); 
hold on 
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('ML method curve fit of %i Decay muon decays in 50 ns bins, third set of 500',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(5),'ML_Curve_Fit_500c.jpg');
%Fourth Set
%Divide the data up into 50 ns bins
maxdecay=max(decaydata4);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(decaydata4,edges);
sz=size(edges);
count=sz(1,2);
bins=zeros(count-1,1);
i=1;
while i<count
    bins(i)=sum(databinned==i);
    i=i+1;
end
binerror=bins.^(1/2);
m=size(bins);
m=m(1,1);
ploterror=zeros(m,1);
j=1;
while j<m
    ploterror(j,1)=binerror(j,1);
    j=j+5;
end
%We assume the decay is somewhere between 0 and 5 microseconds
xL=linspace(0,5,10000)';
n=5/10000;
%Perform the ML to best fit the data with an exponential decay curve
%Determine ln(L) for tau between 0 and 5
lnL=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    expected=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        y(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        i=i+1;
    end
    lnL(j,1)=sum(bins.*log(y*0.05)-(y*0.05));
    j=j+1;
end
%Determine the highest value of ln(L) and set tau to that value
[~,max1index]=min(abs(lnL-(maxlnL-2)));
tau=xL(maxindex,1);
dtau=xL(max1index,1);
error=abs(tau-dtau);
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
figure (6)
set(figure(6),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.'); 
hold on 
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('ML method curve fit of %i Decay muon decays in 50 ns bins, fourth set of 500',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(6),'ML_Curve_Fit_500d.jpg');
%Fifth Set
%Divide the data up into 50 ns bins
maxdecay=max(decaydata5);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(decaydata5,edges);
sz=size(edges);
count=sz(1,2);
bins=zeros(count-1,1);
i=1;
while i<count
    bins(i)=sum(databinned==i);
    i=i+1;
end
binerror=bins.^(1/2);
m=size(bins);
m=m(1,1);
ploterror=zeros(m,1);
j=1;
while j<m
    ploterror(j,1)=binerror(j,1);
    j=j+5;
end
%We assume the decay is somewhere between 0 and 5 microseconds
xL=linspace(0,5,10000)';
n=5/10000;
%Perform the ML to best fit the data with an exponential decay curve
%Determine ln(L) for tau between 0 and 5
lnL=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    expected=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        y(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        i=i+1;
    end
    lnL(j,1)=sum(bins.*log(y*0.05)-(y*0.05));
    j=j+1;
end
%Determine the highest value of ln(L) and set tau to that value
[~,max1index]=min(abs(lnL-(maxlnL-2)));
tau=xL(maxindex,1);
dtau=xL(max1index,1);
error=abs(tau-dtau);
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
figure (7)
set(figure(7),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.'); 
hold on 
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('ML method curve fit of %i Decay muon decays in 50 ns bins, fifth set of 500',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(7),'ML_Curve_Fit_500e.jpg');
%Sixth Set
%Divide the data up into 50 ns bins
maxdecay=max(decaydata6);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(decaydata6,edges);
sz=size(edges);
count=sz(1,2);
bins=zeros(count-1,1);
i=1;
while i<count
    bins(i)=sum(databinned==i);
    i=i+1;
end
binerror=bins.^(1/2);
m=size(bins);
m=m(1,1);
ploterror=zeros(m,1);
j=1;
while j<m
    ploterror(j,1)=binerror(j,1);
    j=j+5;
end
%We assume the decay is somewhere between 0 and 5 microseconds
xL=linspace(0,5,10000)';
n=5/10000;
%Perform the ML to best fit the data with an exponential decay curve
%Determine ln(L) for tau between 0 and 5
lnL=zeros(10000,1);
j=1;
while j<10001
    i=1;
    x=zeros(m,1);
    y=zeros(m,1);
    expected=zeros(m,1);
    while i<m+1
        x(i,1)=0.025+i*.05;
        y(i,1)=(nm*0.05)/((j*n))*exp(-(0.025+i*.05)/(j*n));
        i=i+1;
    end
    lnL(j,1)=sum(bins.*log(y*0.05)-(y*0.05));
    j=j+1;
end
%Determine the highest value of ln(L) and set tau to that value
[max1,max1index]=min(abs(lnL-(maxlnL-2)));
tau=xL(maxindex,1);
dtau=xL(max1index,1);
error=abs(tau-dtau);
%Plot the data and the curve fit
i=1;
x=zeros(m,1);
y=zeros(m,1);
while i<m+1
      x(i,1)=0.025+i*.05;
      y(i,1)=(nm*0.05)/((tau))*exp(-(0.025+i*.05)/(tau));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
figure (8)
set(figure(8),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins, ploterror,'.'); 
hold on 
plot(x,y);
axis([0 20 0 inf]);
caption=sprintf('Curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('ML method curve fit of %i Decay muon decays in 50 ns bins, sixth set of 500',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Decay Data',caption);
saveas(figure(8),'ML_Curve_Fit_500f.jpg');