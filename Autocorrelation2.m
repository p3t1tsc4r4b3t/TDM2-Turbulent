clear all


file_name=['signal.txt'];
data=load(file_name);

n=data(:,1);
t=data(:,2);
U=data(:,3);
V=data(:,4);
W=data(:,5);


T = max(t);
tau = T/10;

%Bruit blanc

bb = randn(length(t),1);
i = xcorr(bb,"unbiased");
i = i(floor(length(i)/2) + 1 : length(i),1);


figure (1)
subplot(2,1,1);plot(t,bb)
xlabel('temps')
ylabel('Amplitude')
grid on

subplot(2,1,2);plot(t,i)
xlabel('temps')
ylabel('auto-corrélation')
grid on

% Signal sinusoidal

f=0.01; % fréquence
C= cos(2*pi.*t*f);

i = xcorr(C,"unbiased");
i = i(floor(length(i)/2)+1:length(i),1);


figure (2)
subplot(2,1,1);plot(t,C)
xlabel('temps')
ylabel('Amplitude')
grid on

subplot(2,1,2);plot(t,i)
xlabel('temps')
ylabel('auto-corrélation')
grid on

% Langevin
for y=3:9
dt=0.0050;
ll= Langevin(0,1,T,dt,length(t));

i = xcorr(ll,"unbiased");
i = i(floor(length(i)/2)+1:length(i),1);


figure (y)
subplot(2,1,1);plot(t,ll)
xlabel('temps')
ylabel('Amplitude')
grid on

subplot(2,1,2);plot(t,i)
xlabel('temps')
ylabel('auto-corrélation')
grid on
end



function X = Langevin(Xmean, Xvar, T, dt, N)
	
	%return a signal given by the Langevin process
	% with:
	% * Xmean: the mean of the process
	% * Xvar: its variance
	% * T:its correlation time 
	% * dt: the time step
	% and N the number of time step
	
	dt_adim=dt/T;
	h=sqrt(Xvar*dt_adim);
	
	X=zeros(N,1);
	X(1)=randn()*sqrt(Xvar);
    for i=2:N
		dx = -(X(i-1) - Xmean) * dt_adim;
		dx = dx + randn()* h;
		X(i) = X(i-1) + dx ;
    end
	
end

