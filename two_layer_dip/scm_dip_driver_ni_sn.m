%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

%Elastic parameters for the top material (nickel)
par.cs = mp('3000'); %transverse wave speed m/sec
par.cl = mp('6040'); %Longitudinal wave Speed m/sec
par.rho = mp('8900'); %Density kg/m^3
par.ell = mp('50*10^-9');
par.c11 = par.rho*par.cl^2;
par.c12 = par.rho*(par.cl^2-2*par.cs^2);
par.c44 = par.rho*par.cs^2;

%Magnetic parameters for the top material
par.h0 = mp('.65');     %applied field
par.ms = mp('480000'); %saturation magnetization
par.exl = mp('7.72')*mp('10^-9');
par.lam = par.exl.^2; %exchange length squared
par.mu0 = 4*mp('pi')*mp('10^-7'); %Vacuum permeability
par.g = mp('1.760859644')*mp('10^11'); %gyromagnetic ratio of an electron
par.b2 = mp('10000000'); %Magnetoelastic coupling constant


%Physical parameters for the bottom material (SiN)
par.CL = mp('8000');
par.CS = mp('3800');
par.RH = mp('3250');
par.C11 = par.RH*par.CL^2;
par.C12 = par.RH*(par.CL^2-2*par.CS^2);
par.C44 = par.RH*par.CS^2;
par.ELL = mp('100*10^-9');

N = 12;

%Chebyshev differentiation matrices
[D,~] = cheb(N-1);
Dt = 2*D./par.ell;
Db = 2*D./par.ELL;
Dt2 = Dt*Dt;
Db2 = Db*Db;

%Iterate omega and k
% k = om/par.cl;
n = 360;
kmin = mp('3*10^8');
kmax = mp('4.2*10^8');
%kmax = 2*mp('pi')*wmax/(par.cs*mp('.9'));
%kmax = mp('.35')*mp('10^8');

%kv = omv./par.cs;
kv = linspace(kmin,kmax,n);

%kv = mp('45 *10^7');

% omv =     mp('158263157894.7368421052631578947368');
% kv = mp('.1');

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

    wlong_top = wm;
    wshear_top = wm;
    wlong_bot = wm;
    wshear_bot = wm;
    wmag = wm;
    
    wlong_top(tp~=1)=NaN;
    wshear_top(tp~=2)=NaN;
    wlong_bot(tp~=3) = NaN;
    wshear_bot(tp~=4) = NaN;
    wmag(tp~=5)=NaN;
   
    kplot = double(kv);
    
      data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\double_layer_dip_ni_sn_N_12.mat'];
     save(data_dir,'par','wlong_top','wlong_bot','wshear_top','wshear_bot','wmag','kplot','N');

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
    plot(kplot/10^9,sort(real(wshear_top))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wshear_bot))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');


    hold off
    % labels
    axis xy
    xlabel('Wavenumber [nm$^{-1}$]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');


 axis([3*10^8 4*10^8 120*10^9 240*10^9]/10^9)
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
set(gcf,'Resize','off')


%print(fh,'-dpdf','two_layer_dip_ni_sn_N_6_v2');
