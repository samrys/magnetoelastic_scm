%%Runs SCM solver for magnetoelastic waves in one layer
%   Includes dipole field
%   Calls cheb_mat_dip_one_layer.m and mag_class_one_layer.m

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

%Elastic parameters for the top material (currently YIG)
par.cs = mp('3800'); %transverse wave speed m/sec
par.cl = mp('7200'); %Longitudinal wave Speed m/sec
par.rho = mp('5170'); %Density kg/m^3
par.ell = mp('30*10^-9'); %Width of material layer (paper calls this d)
par.c11 = par.rho*par.cl^2;
par.c12 = par.rho*(par.cl^2-2*par.cs^2);
par.c44 = par.rho*par.cs^2;

%Magnetic parameters for the top material
par.ms = mp('140000'); %saturation magnetization
par.lam = mp('3')*mp('10^-16'); %exchange length squared
par.mu0 = 4*mp('pi')*mp('10^-7'); %Vacuum permeability
par.h0 = mp('.25');     %applied field
par.g = mp('1.760859644')*mp('10^11'); %gyromagnetic ratio of an electron
par.b2 = mp('1')*mp('7000000'); %Magnetoelastic coupling constant

N = 16;

%Chebyshev differentiation matrices
[D,~] = cheb(N-1);
Dt = 2*D./par.ell;
Dt2 = Dt*Dt;

%Iterate k
n = 100; %Choose discretization of k
kmin = mp('3.1*10^8');
kmax = mp('4*10^8');
kv = linspace(kmin,kmax,n);

wm = zeros(12*N,length(kv));
tp = wm;

h = @(kvar) cheb_mat_dip_one_layer(kvar,Dt,Dt2,N,par);

parfor ii = 1:length(kv)
    k = kv(ii);
    [A,B,C] = h(k);
    [Xv,eigval] = polyeig(A,B,C);% eigval = polyeig(A,B,C);
    tp(:,ii) = mag_class_one_layer(Xv,eigval,k,N,par);
    wm(:,ii) = eigval;
end
    wlong = wm;
    wtran = wm;
    wmag = wm;
    
    wlong(tp~=1)=NaN;
    wtran(tp~=3)=NaN;
    wmag(tp~=2)=NaN;

    

  kplot = double(kv);
    
%      data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\one_layer_dip_yig_N_16_b2_zero_v2.mat'];
%     save(data_dir,'par','wlong','wmag','wtran','kplot','N');

%      data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data\one_layer_dip_yig_N_12_b2_p1.mat'];
 
 
     %% Create figure
 %    load(data_dir,'par','wlong','wmag','wtran','kplot','N');
     
    %Figure specs
ml = 0.11; % Margin left
mr = 0.03; % Margin right
mt = 0.04; % Margin top
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
    
   
    hold on


    plot(kplot/10^9,sort(real(wtran))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');

    hold off
    % labels
    axis xy
    xlabel('Wavenumber [rad/nm]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');

 axis([3.1*10^8 3.8*10^8 170*10^9 250*10^9]/10^9)    
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'linewidth',lw);     
% set(gca,'xtick',[.32 .34 .36 .38],'ytick',[200 220 240])   
set(gcf,'Resize','off')


%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw/2);
plot(nan,nan,'-g','LineWidth',lw/2);
hold off

legend('Shear','Magnetic',...
    'interpreter','latex','location','SE','fontsize',fontsize)
legend boxoff 
set(gcf,'Resize','on')


%print(fh,'-dpdf','one_layer_dip_yig_30_b2_nonzero_v2_new_2.pdf');
