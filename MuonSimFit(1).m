%Get the number of decays to be simulated
nm=inputdlg('How many muons decays do you want to simulate?', 'Simulation Size');
nm=string(nm);
nm=str2double(nm);
%Simulate the muon decays
i=1;
simdata=zeros(nm,1);
while i<nm+1
    simdata(i,1)=-2.2*log(rand);
    i=i+1;
end
%Binned Maximum Likelyhood method fitting
%Divide the data up into 50 ns bins
maxdecay=max(simdata);
lastedge=ceil(maxdecay);
edges=0:.05:lastedge;
databinned=discretize(simdata,edges);
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
[max1,max1index]=min(abs(lnL-(maxlnL-2)));
tauML=xL(maxindex,1);
dtauML=xL(max1index,1);
errorML=abs(dtauML-tauML);
%Plot the data and the curve fit
i=1;
xML=zeros(m,1);
yML=zeros(m,1);
while i<m+1
      xML(i,1)=0.025+i*.05;
      yML(i,1)=(nm*0.05)/((tauML))*exp(-(0.025+i*.05)/(tauML));
      i=i+1;
end
edges=transpose(edges);
edges(end)=[];
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
figure (1)
set(figure(1),'Units','centimeters','Position',[0 0 25 15]);
errorbar(edges,bins,ploterror,'.');
hold on
plot(x,y, xML, yML);
axis([0 20 0 inf]);
captionML=sprintf('ML curve fit with mean decay time (microseconds) %.2f +/- %.2f',tauML,errorML);
caption=sprintf('LS curve fit with mean decay time (microseconds) %.2f +/- %.2f',tau,error);
tcaption=sprintf('Curve fits of %i simulated muon decays in 50 ns bins.',nm);
title(tcaption);
xlabel('Decay time in us');
ylabel('Number of decays');
legend('Simulated Data',caption,captionML);



