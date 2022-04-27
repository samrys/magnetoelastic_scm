%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

%Elastic parameters for the top material
par.cs = mp('3800'); %transverse wave speed m/sec
par.cl = mp('7200'); %Longitudinal wave Speed m/sec
par.rho = mp('5170'); %Density kg/m^3
par.ell = mp('30*10^-9');
par.c11 = par.rho*par.cl^2;
par.c12 = par.rho*(par.cl^2-2*par.cs^2);
par.c44 = par.rho*par.cs^2;

%Magnetic parameters for the top material
par.h0 = mp('.25');     %applied field
par.ms = mp('140000'); %saturation magnetization
par.lam = mp('3')*mp('10^-16'); %exchange length squared
par.mu0 = 4*mp('pi')*mp('10^-7'); %Vacuum permeability
par.g = mp('1.760859644')*mp('10^11'); %gyromagnetic ratio of an electron
par.b2 = mp('7000000'); %Magnetoelastic coupling constant


%Physical parameters for the bottom material
par.CL = mp('6400');
par.CS = mp('3500');
par.RH = mp('7080');
par.C11 = par.RH*par.CL^2;
par.C12 = par.RH*(par.CL^2-2*par.CS^2);
par.C44 = par.RH*par.CS^2;
par.ELL = mp('50*10^-9');

N = 8;

%Chebyshev differentiation matrices
[D,~] = cheb(N-1);
Dt = 2*D./par.ell;
Db = 2*D./par.ELL;
Dt2 = Dt*Dt;
Db2 = Db*Db;

%Iterate omega and k
% k = om/par.cl;
n = 50;
kmin = mp('3*10^8');
kmax = mp('4*10^8');
kv = linspace(kmin,kmax,n);



wm = zeros(18*N,length(kv));
tp = wm;
st = wm;

h = @(kvar) cheb_mat_dip(kvar,Dt,Db,Dt2,Db2,N,par);

parfor ii = 1:length(kv)
    k = kv(ii);
    [A,B,C] = h(k);
    [Xv,eigval] = polyeig(A,B,C);   
    tp(:,ii) = mag_class(Xv,eigval,k,N,par);
    wm(:,ii) = eigval;
end

    wtop1 = wm;
    wtop2 = wm;
    wtop3 = wm;
    wbot1 = wm;
    wbot2 = wm;
    wbot3 = wm;
    wmag = wm;
    
    wtop1(tp~=1) = NaN;
    wtop2(tp~=2) = NaN;
    wtop3(tp~=3) = NaN;
    wbot1(tp~=4) = NaN;
    wbot2(tp~=5) = NaN;
    wbot3(tp~=6) = NaN;
    wmag(tp~=7)  = NaN;
   
    kplot = double(kv);
    
%     data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\double_layer_dip_yig_ggg_N_12_zoom.mat'];
%    save(data_dir,'par','wtop1','wtop2','wtop3','wbot1','wbot2','wbot3','wmag','kplot','N');

 %% plot results 

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
    plot(kplot/10^9,sort(real(wtop2))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
        plot(kplot/10^9,sort(real(wtop3))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wbot2))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
        plot(kplot/10^9,sort(real(wbot3))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');


    hold off
    
    % labels
    axis xy
    xlabel('Wavenumber [nm$^{-1}$]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');


 axis([3.2*10^8 4*10^8 180*10^9 250*10^9]/10^9)
%     set(gca,'ytick',[230,250,270,290],...
%         'xtick',[0,.36,.37,.38,.39,.4,.41,.42,.43,.44,.45,.6]);  
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex');     
set(gcf,'Resize','off')


%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw);
plot(nan,nan,'-r','LineWidth',lw);
plot(nan,nan,'-g','LineWidth',lw);
hold off

legend('Top Shear','Bottom Shear','Magnetic',...
    'interpreter','latex','location','SE','fontsize',fontsize)
legend boxoff 
 set(gcf,'Resize','on')


%print(fh,'-dpdf','two_layer_dip_yig_sn_N_6_v2');
