T = readtable('USGS2015-2020.csv');
discharge = T.x293221_00060_00003;
time = T.datetime;


%hist(discharge)
%mean(discharge)

n = 200; %length(time);
dt = 1;
t = 0:dt:n-1;
a = 1.32; %flow decay
b = 1.53; %brownian turbulance
sr = 0.44; %saturation decay rate
sp = 0.37; %saturation percentage
lambda = 0.86; %precipitation frequency
M = 340; %precipitation magnitude coefficient

dBt = normrnd(0,1,[1,n])
dPt = poissrnd(lambda)

X = zeros(1,n);
s = zeros(1,n); s(1) = 100;
dx = 0; ds = 0;
for j = 2:n
    N = poissrnd(lambda)
    V = exprnd(M,[1,N])
    dx = (s(j-1)-a*X(j-1))*dt + b*dBt(j) + sum(V);
    ds = -(sr*s(j-1))*dt + sp*sqrt(X(j-1))*dt
    X(j) = X(j-1) + dx;
    s(j) = s(j-1) + ds;
end

%tiledlayout(2,1)
%nexttile
plot(t,X);%hold on
%plot(t,s,"--");hold off
%nexttile

%plot(t,X)
xlabel('Time') 
ylabel('Discharge (ft^3 / s)') 
title('Discharge model: With Saturation');
mean(X)
mean(discharge)
