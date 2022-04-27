%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

%Elastic parameters for the top material
par.cs = mp('2900'); %transverse wave speed m/sec
par.cl = mp('5500'); %Longitudinal wave Speed m/sec
par.rho = mp('8900'); %Density kg/m^3
par.ell = mp('30*10^-9');
par.c11 = par.rho*par.cl^2;
par.c12 = par.rho*(par.cl^2-2*par.cs^2);
par.c44 = par.rho*par.cs^2;

%Magnetic parameters for the top material
%par.h0 = mp('.5');     %applied field
par.h0 = mp('.2');     %applied field
par.ms = mp('480000'); %saturation magnetization
par.lam = mp('5.929')*mp('10^-17'); %exchange length squared
par.mu0 = 4*mp('pi')*mp('10^-7'); %Vacuum permeability
par.g = mp('1.760859644')*mp('10^11'); %gyromagnetic ratio of an electron
par.b2 =mp('10000000'); %Magnetoelastic coupling constant
%par.b2 = mp('0');
par.beta = par.g.*par.mu0.*par.lam.*par.ms;
par.eps = 1/par.cs*sqrt(par.g*par.beta*par.b2.^2/(par.c44*par.ms));

N = 6;

%Chebyshev differentiation matrices
[D,~] = cheb(N-1);
Dt = 2*D./par.ell;
Dt2 = Dt*Dt;

%Iterate omega and k
% k = om/par.cl;
n = 300;
kmin = mp('4*10^8');
kmax = mp('5*10^8');
%kmax = 2*mp('pi')*wmax/(par.cs*mp('.9'));
%kmax = mp('.35')*mp('10^8');

%kv = omv./par.cs;
kv = linspace(kmin,kmax,n);

%kv = mp('45 *10^7');

% omv =     mp('158263157894.7368421052631578947368');
% kv = mp('.1');

wm = zeros(6*N,length(kv));
tp = wm;

h = @(kvar) cheb_mat_one_layer_simp(kvar,Dt,Dt2,N,par);
par2 = par;
par2.h0 = mp('1');
h2 = @(kvar) cheb_mat_one_layer_simp(kvar,Dt,Dt2,N,par2);

for ii = 1:length(kv)
    k = kv(ii);
    [A,B,C] = h(k);
    [Xv,eigval] = polyeig(A,B,C);% eigval = polyeig(A,B,C);
    tp(:,ii) = mag_class_one_layer_simp(Xv,eigval,k,N,par);
    wm(:,ii) = eigval;
end
%kplot_mag = double(kv);
%wplot_mag = double(wm);
%save('C:\Users\samry\Documents\MATLAB\Elastic\data\fixed_h0_N_12_one_layer_d_5nm.mat','kplot_mag','wplot_mag','par');
%    data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\data\stor_',num2str(n),'_no_mag.mat'];
    %save(data_dir,'par','wmin','wmax','kmin','kmax','A','n');

 %% plot results 

   %Figure specs
ml = 0.11; % Margin left
mr = 0.03; % Margin right
mt = 0.03; % Margin top
mb = 0.17;  % Margin bottom
pb = 0.08; % Interaxes padding bottom
pr = 0.085; % Interaxes padding right
Nplots = 1;
spanx = (1-ml-mr-(Nplots-1)*pr)/Nplots;
spany = (1-mt-mb-(1-1)*pb)/1;
fac = 4; % Factor to increase figure size (dashed line hack)
fig_width = 8.6*fac; % in cm
fig_height = 4.5*fac;
fontsize = 8*fac;
lw = 1.1*fac;

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

set(gcf,'Resize','off')

axes('Position',[ml,mb,spanx,spany]);
    
    
    wlong = wm;
    wtran = wm;
    wmag = wm;
    
    wlong(tp~=1)=NaN;
    wtran(tp~=3)=NaN;
    wmag(tp~=2)=NaN;

    
    kplot = double(kv);
    hold on

    plot(kplot/10^9,sort(real(wlong))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wtran))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');

    hold off
    % labels
    axis xy
    xlabel('Wavenumber [rad/nm]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex','xtick',[.42 .44 .46 .48 .5]);


axis([.42*10^9 .5*10^9 190*10^9 260*10^9]/10^9)
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'linewidth',fac/2);     
set(gcf,'Resize','off')


nl            =   5500        ;
ns            =   2900        ;

wnl = nl*kplot/(2*pi);
wns = ns*kplot/(2*pi);
wmag = par.g.*(par.h0+par.mu0*par.lam*(kplot).^2*par.ms)/(2*pi);

hold on
% plot(kplot/10^9,wnl/10^9,'--k','LineWidth',lw,'HandleVisibility','off');
% plot(kplot/10^9,wns/10^9,'--k','LineWidth',lw,'HandleVisibility','off');

% plot(kplot/10^9,wmag/10^9,'--k','LineWidth',lw,'HandleVisibility','off');
hold off

%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw);
plot(nan,nan,'-g','LineWidth',lw);
plot(nan,nan,'r-.','Linewidth',lw)

hold off

legend('Shear','Magnetic','Asymptotic',...
    'interpreter','latex','location','SE','fontsize',fontsize)
legend boxoff 

set(gcf,'Resize','off')


K = par.cs/par.beta;
Om = par.cs*K;
H = par.g*par.h0*par.beta/par.cs^2;
G = par.cl/par.cs;
%Plot asymptotic prediction
for nn = 0:1
an = nn*par.beta*pi/(par.cs*par.ell);
kstar = sqrt(-H+1/2*(1-2*an^2+sqrt(1-4*H+4*an^2*(G^2-1))));
wstar = sqrt(kstar^2+G^2*an^2);
del = (kplot/K-kstar)/par.eps;
kasymp = K*(kstar+par.eps*del);
fasymp1 = Om*(wstar+par.eps.*kstar/2*(del.*(1/(wstar)+2)+...
     sqrt(del.^2*(2*wstar-1)^2+wstar)))/(2*pi);
fasymp2 = Om*(wstar+par.eps.*kstar/2*(del.*(1/(wstar)+2)-...
     sqrt(del.^2*(2*wstar-1)^2+wstar)))/(2*pi);
%fasymp1 = Om/2*(3*del+sqrt(del.^2+1))/(2*pi);
%fasymp2 = Om/2*(3*del-sqrt(del.^2+1))/(2*pi);

hold on
plot(kplot/10^9,fasymp1/10^9,'r-.','Linewidth',1.5*lw,'HandleVisibility','off')
plot(kplot/10^9,fasymp2/10^9,'r-.','Linewidth',1.5*lw,'HandleVisibility','off')
% if nn>-1
% mag1 = Om*(kplot.^2/K^2+H+an^2);
% el1 = Om*sqrt(kplot.^2/K^2+G^2*an^2);
% plot(kplot/10^9,mag1/(2*pi*10^9),'r-.','Linewidth',lw)
% plot(kplot/10^9,el1/(2*pi*10^9),'r-.','Linewidth',lw)
% end
plot(kstar*K/10^9,wstar*Om/(2*pi*10^9),'*k','markersize',4*lw,'handlevisibility','off')

hold off
end
    text(.4485,205,'$(n,j)=(0,0)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',8*fac)
    text(.463,240,'$(n,j)=(1,1)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','right','fontsize',8*fac)
print(fh,'-dpdf','asymptotic_comparison.pdf');
