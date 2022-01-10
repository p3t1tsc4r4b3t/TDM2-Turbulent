close all 
signal = load('signal.txt'); 

t=signal(:,2) ;
U=[signal(:,3) signal(:,4) signal(:,5)] ;

% 1.courbe de u(t)  
%plot(t,U);      

% 2.calcul de la moyenne temporelle U_moy
% U_moy : Tableau à 3 elements (moyenne selon x, y et z)
U_moy = zeros(1,3) ;
U_moy(1) = mean(U(:,1)) ;
U_moy(2) = mean(U(:,2)) ;
U_moy(3) = mean(U(:,3)) ;


% 2.calcul de l'ecart-type U_std
% U_std : Tableau à 3 elements (ecart-type selon x, y et z)
% on a bien pris ici la vitesse fluctuante (compris dans std)
U_std = zeros(1,3) ;
U_std(1) = std(U(:,1)) ;
U_std(2) = std(U(:,2)) ;
U_std(3) = std(U(:,3)) ;

%3.Calcul de la densite de probabilite de U: U_prob
U_prob = zeros(100,3);
U_i = zeros(100,3);
[U_prob(:,1),U_i(:,1)]=ksdensity(U(:,1));
[U_prob(:,2),U_i(:,2)]=ksdensity(U(:,2));
[U_prob(:,3),U_i(:,3)]=ksdensity(U(:,3));
%figure ; 

%plot(U_i(:,1),U_prob(:,1)) ;

%3.On verifie que la somme de la probabilité vaut bien 1
%Somme_proba=[trapz(U_i(:,1),U_prob(:,1)) trapz(U_i(:,2),U_prob(:,2)) trapz(U_i(:,3),U_prob(:,3))];

%3.Calcul de la loi normale
Gauss_prob = zeros(100,3) ;
Gauss_prob(:,1)= normpdf(U_i(:,1),U_moy(:,1),U_std(:,1));
Gauss_prob(:,2)= normpdf(U_i(:,2),U_moy(:,2),U_std(:,2));
Gauss_prob(:,3)= normpdf(U_i(:,3),U_moy(:,3),U_std(:,3));

% On trace les courbes que l'on compare \ufffd une loi normale N(U_moy, U_std)
%figure;
%plot(U_i(:,1),U_prob(:,1));
% hold on;
% plot(U_i(:,1),Gauss_prob(:,1));

% % 3.Comparaison des courbes de densite (avec les gaussiennes) pour u,v et w
% plot(U_i(:,1),U_prob(:,1), U_i(:,1),Gauss_prob(:,1)) ;
% legend('U', 'U Gaussien')
% xlabel('Vitesse')
% ylabel('Densite de probabilite')
% figure;
% plot(U_i(:,2),U_prob(:,2), U_i(:,2),Gauss_prob(:,2)) ;
% legend('V', 'V Gaussien')
% xlabel('Vitesse')
% ylabel('Densite de probabilite')
% figure;
% plot(U_i(:,3),U_prob(:,3), U_i(:,3),Gauss_prob(:,3)) ;
% legend('W', 'W Gaussien')
% xlabel('Vitesse')
% ylabel('Densite de probabilite')

% 3.Second calcul de la moyenne et de l'ecart-type grace a la densite de proba
U_moy2 = [trapz(U_i(:,1),U_prob(:,1).*U_i(:,1)) trapz(U_i(:,2),U_prob(:,2).*U_i(:,2)) trapz(U_i(:,3),U_prob(:,3).*U_i(:,3))];
U_std2 = [sqrt(trapz(U_i(:,1),U_prob(:,1).*(U_moy2(1) -U_i(:,1)).*(U_moy2(1) -U_i(:,1)))) sqrt(trapz(U_i(:,2),U_prob(:,2).*(U_moy2(2) -U_i(:,2)).*(U_moy2(2) -U_i(:,2)))) sqrt(trapz(U_i(:,3),U_prob(:,3).*(U_moy2(3) -U_i(:,3)).*(U_moy2(3) -U_i(:,3))))];
% On a bien U_moy2 ~ U_moy et U_std2 ~ U_std 

% 4. Calcul de la d\ufffdriv\ufffde temporelle du signal de vitesse
U_der = zeros(length(U(:,1))-1,3);
for i=1:length(U_der(:,1))
   U_der(i,1) = (U(i+1,1)-U(i,1))/(t(i+1)-t(i));
   U_der(i,2) = (U(i+1,2)-U(i,2))/(t(i+1)-t(i));
   U_der(i,3) = (U(i+1,3)-U(i,3))/(t(i+1)-t(i));
end
% precautions: On perd une cellule de calcul car il faut N+1 donn\ufffdes pour
% calculer N d\ufffdriv\ufffdes


% 5. Calcul de la moyenne/densit\ufffd de cette d\ufffdriv\ufffde
% U_der_moy : Tableau a 3 elements (moyenne selon x, y et z)
U_der_moy = zeros(1,3) ;
U_der_moy(1) = mean(U_der(:,1)) ;
U_der_moy(2) = mean(U_der(:,2)) ;
U_der_moy(3) = mean(U_der(:,3)) ;
    
% calcul de l'ecart-type U_std
% U_der_std : Tableau a 3 elements (ecart-type selon x, y et z)
U_der_std = zeros(1,3) ;
U_der_std(1) = std(U_der(:,1)) ;
U_der_std(2) = std(U_der(:,2)) ;
U_der_std(3) = std(U_der(:,3)) ;

% Calcul de la densit\ufffd de probabilit\ufffd de la d\ufffdriv\ufffde temporelle de la
% vitesse
U_der_prob = zeros(100,3);
U_der_i = zeros(100,3);
[U_der_prob(:,1),U_der_i(:,1)]=ksdensity(U_der(:,1));
[U_der_prob(:,2),U_der_i(:,2)]=ksdensity(U_der(:,2));
[U_der_prob(:,3),U_der_i(:,3)]=ksdensity(U_der(:,3));

% Calcul de la loi normale
Gauss_der_prob = zeros(100,3) ;
Gauss_der_prob(:,1)= normpdf(U_der_i(:,1),U_der_moy(:,1),U_der_std(:,1));
Gauss_der_prob(:,2)= normpdf(U_der_i(:,2),U_der_moy(:,2),U_der_std(:,2));
Gauss_der_prob(:,3)= normpdf(U_der_i(:,3),U_der_moy(:,3),U_std(:,3));

% figure;
% plot(U_der_i(:,1),U_der_prob(:,1));
% hold on;
% plot(U_der_i(:,1),Gauss_der_prob(:,1));

% Comparaison des densit\ufffds de probabilit\ufffd U vs dU/dt:
% plot(U_der_i(:,1),U_der_prob(:,1), 'R');
% hold on;
% plot(U_i(:,1),U_prob(:,1), 'G');

% 
% % 2.3 Autocorrelation temporelle du signal
% %1.Signal bruit blanc
% s_bc_0= zeros(1000,1);
% s_bc_10= zeros(1000,1);
% for k=1:1000
%    s_bc_0(k)= 2*rand()-1 ;
%    s_bc_10(k)= 10*rand() ;
% end
%  [rhobc,lagsbc] = xcorr(s_bc_0(:,1), "normalized");
%  [rhobc_10,lagsbc_10] = xcorr(s_bc_10(:,1), "normalized");
%  rhobc = rhobc(round(size(rhobc)/2):size(rhobc));
%  a = size(lagsbc);
%  lagsbc = lagsbc(round(a(2)/2):a(2));
%  lagscbc = lagsbc(1:round(a(2)/20));
%  rhocbc =  rhobc(1:round(a(2)/20));
%  
% 
%  rhobc_10 = rhobc_10(round(size(rhobc_10)/2):size(rhobc_10));
%  a = size(lagsbc_10);
%  lagsbc_10 = lagsbc_10(round(a(2)/2):a(2));
%  lagscbc_10 = lagsbc_10(1:round(a(2)/20));
%  rhocbc_10 =  rhobc_10(1:round(a(2)/20));
%  
%  plot(lagscbc, rhocbc,lagscbc_10, rhocbc_10);
% %On prend seulement les valeurs de tau < T/10
% %Calcul de l'int\ufffdgale Tau int
%  tau_intbc = trapz(lagscbc, rhocbc) ;
%  
% figure;
% 
% % 2.3 Autocorrelation temporelle du signal
% %1.Signal sinusoidal
% s_sin= zeros(1000,1);
% for k=1:1000
%    s_sin(k)= sin(pi*k/8);
% end
%  [rhosin,lagssin] = xcorr(s_sin(:,1), "normalized");
%  rhosin = rhosin(round(size(rhosin)/2):size(rhosin));
%  a = size(lagssin);
%  lagssin = lagssin(round(a(2)/2):a(2));
%  lagscsin = lagssin(1:round(a(2)/20));
%  rhocsin =  rhosin(1:round(a(2)/20));
%  plot(lagscsin, rhocsin);
% %On prend seulement les valeurs de tau < T/10
% %Calcul de l'int\ufffdgale Tau int
%  tau_intsin = trapz(lagscsin, rhocsin) ;
%  


%Xmean = 1
%Xvar1 = 1
%Xvar2 = 100
%Xvar3 = 1000
%T = 0.1
	
%s_lg = zeros(1000,1);
%s_lg= Langevin(Xmean,Xvar1,T,0.02,1000);
%s_lg2 = zeros(1000,1);
%s_lg2 = Langevin(Xmean,Xvar2,T,0.02,1000);
%s_lg3 = zeros(1000,1);
%s_lg3 = Langevin(Xmean,Xvar3,T,0.02,1000);


%Xmean1 = 1
%Xmean2 = 10
%Xmean3 = 100
%Xvar = 10
%T = 0.1
	
%s_lg = zeros(1000,1);
%s_lg= Langevin(Xmean1,Xvar,T,0.02,1000);
%s_lg2 = zeros(1000,1);
%s_lg2 = Langevin(Xmean2,Xvar,T,0.02,1000);
%s_lg3 = zeros(1000,1);
%s_lg3 = Langevin(Xmean3,Xvar,T,0.02,1000);

% 2.3 Autocorrelation temporelle du signal
%1.Signal de Langevin
%[rholg1,lagslg1] = xcorr(s_lg(:,1), "unbiased");
%[rholg2,lagslg2] = xcorr(s_lg2(:,1), "unbiased");
%[rholg3,lagslg3] = xcorr(s_lg3(:,1), "unbiased");
%tau = 100;
%i = 0;
%j=0;
%rho = zeros(tau,1);
%lags = linspace(0,tau,tau);
%for i=1:tau
%    for j=1:(length(s_lg)-i)
%        rho(i,1)=rho(i,1)+(s_lg(j)-mean(s_lg))*(s_lg(j+i-1)-mean(s_lg));
%    end
%rho(i,1)=rho(i,1)/(length(s_lg)-i-1);
%end
%rho(:,1)=rho(:,1)/mean((s_lg-mean(s_lg)).^2);
%lagsa = lags;
%rhoa = rho;

 % rholg = rholg(round(size(rholg)/2):size(rholg));
  %a = size(lagslg);
  %lagslg = lagslg(round(a(2)/2):a(2));
  %lagsclg = lagslg(1:round(a(2)/20));
  %rhocg =  rhocg(1:round(a(2)/20));%

%rho2 = zeros(tau,1);
%lags2 = linspace(0,tau,tau);
%for i=1:tau
%    for j=1:(length(s_lg2)-i)
%        rho2(i,1)=rho2(i,1)+(s_lg2(j)-mean(s_lg2))*(s_lg2(j+i-1)-mean(s_lg2));
%    end
%rho2(i,1)=rho2(i,1)/(length(s_lg2)-i-1);
%end
%rho2(:,1)=rho2(:,1)/mean((s_lg2-mean(s_lg2)).^2);
%lagsa2 = lags2;
%rhoa2 = rho2;

 % rholg = rholg(round(size(rholg)/2):size(rholg));
 % a = size(lagslg);
  %lagslg = lagslg(round(a(2)/2):a(2));
 % lagsclg = lagslg(1:round(a(2)/20));
 % rhocg =  rhocg(1:round(a(2)/20));%

%rho3 = zeros(tau,1);
%lags3 = linspace(0,tau,tau);
%for i=1:tau
%    for j=1:(length(s_lg3)-i)
%       rho3(i,1)=rho3(i,1)+(s_lg3(j)-mean(s_lg3))*(s_lg3(j+i-1)-mean(s_lg3));
%    end
%rho3(i,1)=rho3(i,1)/(length(s_lg3)-i-1);
%end
%rho3(:,1)=rho3(:,1)/mean((s_lg3-mean(s_lg3)).^2);
%lagsa3 = lags3;
%rhoa3 = rho3;
%
%  rholg = rholg(round(size(rholg)/2):size(rholg));
%  a = size(lagslg);
%  lagslg = lagslg(round(a(2)/2):a(2));
%  lagsclg = lagslg(1:round(a(2)/20));
%  rhocg =  rhocg(1:round(a(2)/20));%


 %plot(lags,rho,lags2,rho2,lags3,rho3);
%xlabel('tau(s)')
%ylabel('autocorrection')
%title('autocorrection avec la variance de Xvari')


% lags100 = lags;
% rho100 = rho;
% On prend seulement les valeurs de tau < T/10
%Calcul de l'int\ufffdgale Tau int
% tau_intlg = trapz(lagsclg, rhocg) ;

  %figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODIFIER LES VALEURS DE LA FCT LANGEVIN POUR OBTENIR UNE COURBE CARACTERISTIQUE ET EXPLOITABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.Calcul du coeff d'autocorrelation et du tps int\ufffdgral pour la vitesse U
%(meme calcul pour V et W)
[rho,lags] = xcorr(U(:,3), "unbiased");

 rho = rho(round(size(rho)/2):size(rho));
 a = size(lags);
 lags = lags(round(a(2)/2):a(2));
 lagsc = lags(1:round(a(2)/20))*0.005;
 rhoc =  rho(1:round(a(2)/20));

 plot(lagsc,rhoc);
 xlabel('W')
 ylabel('autocorrélation')
 grid
 %On prend seulement les valeurs de tau < T/10
%calcul de l'int\ufffdgale Tau int
 tau_int = trapz(lagsc, rhoc) ;
S=tau_int
 % Auto-correlation, courbe non coup\ufffde (pour faire des plots et comparer \ufffd
 % l'hypoth\ufffdse de taylors):
[rho_full,lags_full] = xcorr(U(1:4,1));


%figure;
%plot(lags_full(1,:), rho_full(:,1));
% HYPOTHESE DE TAYLORS
% On suppose que la turb au temps t est la m\ufffdme que la turbulence d\ufffdplac\ufffde
% de U_m/t. On retrouve alors U(x0,t) = U(x1,t+delta_x/U_m).
% On pose dt = dx/U_moy => changement de variable dans l'integrale
% d'autocor.
% On trouve \rho(\tau) = \frac{1}{U_mT-\tau} \int{0}{L__x}\frac
%{u_i^'(x)u_i^'(x+X)dx}{u_i^'^2}

% Donn\ufffdes du fichier:
%# Autocorrelation spatiale dans la direction longitudinale
%# pour les trois composatne de la vitesse
%# rho_ux : autocorrelation pour la composante longitudinale
%# rho_uy : autocorrelation pour la composante normale \ufffd la paroi
%# rho_uz : autocorrelation pour la composante transversale
%# rx : longueur du d\ufffdcalage spatiale
%# ix : nombre de maille du d\ufffdcalage
%# Nombre de stat = 573
%# ix rx rho_ux rho_uy  rho_uz
%#plan i=48
 

%autocorr = load('autocorrection.txt'); 

%ix = autocorr(:,1); 
%rx = autocorr(:,2); 
%rho_x = autocorr(:,3);
%rho_y = autocorr(:,4);
%rho_z = autocorr(:,5);


 %figure;
 % plot(ix, rho_x,ix,rho_y,ix,rho_z);
 
% figure;
%  plot(lags_full, rho_full);
%Signal de lANGEVIN
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
	end;
	
end