%Generate plots for given frequencies and wavenumbers

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

xpts = [.204 .2193 .22573 .24342 .26271 .27719 .29166 .30613 .3206 .3399 .35437 .37367 .38815 .40583 .42352 .44925 .47658];
ypts = [180.4 191.732 198.7907 213.9487 231.2037 243.7736 256.3336 268.9018 281.5757 298.3496 310.9298 327.4389 337.9144 349.4933 358.7379 373.3658 388.452];


kvec    = [mp('.30613*10^9') mp('.22573*10^9')  mp('.44925*10^9')];
fvec = [mp('268.9018*10^9') mp('198.7907*10^9')  mp('373.3658*10^9')];
N = 12;

%Elastic parameters for the top material
par.cs = mp('4000'); %transverse wave speed m/sec
par.cl = mp('5500'); %Longitudinal wave Speed m/sec
par.rho = mp('8900'); %Density kg/m^3
par.ell = mp('50*10^-9');
par.c11 = par.rho*par.cl^2;
par.c12 = par.rho*(par.cl^2-2*par.cs^2);
par.c44 = par.rho*par.cs^2;

%Magnetic parameters for the top material
par.h0 = mp('.25');     %applied field
par.ms = mp('480000'); %saturation magnetization
par.exl = mp('7.72')*mp('10^-9');
par.lam = par.exl.^2; %exchange length squared
par.mu0 = 4*mp('pi')*mp('10^-7'); %Vacuum permeability
par.g = mp('1.760859644')*mp('10^11'); %gyromagnetic ratio of an electron
par.b2 = mp('10000000'); %Magnetoelastic coupling constant


%Physical parameters for the bottom material
par.CL = mp('7700');
par.CS = mp('4600');
par.RH = mp('3170');
par.C11 = par.RH*par.CL^2;
par.C12 = par.RH*(par.CL^2-2*par.CS^2);
par.C44 = par.RH*par.CS^2;
par.ELL = mp('100*10^-9');

ell = double(par.ell)*10^9;
ELL = double(par.ELL)*10^9;

xtop = linspace(0,ell,150);
xbot = linspace(-ELL,0,300);

[s1,fnew1] = wave_structure(fvec(1),kvec(1),par,N,xtop,xbot);
[s2,fnew2] = wave_structure(fvec(2),kvec(2),par,N,xtop,xbot);
[s3,fnew3] = wave_structure(fvec(3),kvec(3),par,N,xtop,xbot);


%str = {'$m$','$A_1$','$A_2$','$A_3$','$\overline{A_1}$','$\overline{A_2}$','$\overline{A_3}$'};

% locvec = [4,4,4];
% txtvec = {'$(a)$','$(b)$','$(c)$'};

zx = -1.05:.05:1.05; %x axis
zy = -105:.5:105;      % y axis

%% Plotting
    %Figure specs
ml = 0.11; % Margin left
mr = 0.02; % Margin right
mt = 0.04; % Margin top
mb = 0.17;  % Margin bottom
pb = 0.08; % Interaxes padding bottom
pr = 0.04; % Interaxes padding right
Nplots = length(fvec);
spanx = (1-ml-mr-(Nplots-1)*pr)/Nplots;
spany = (1-mt-mb-(1-1)*pb)/1;
fac = 4; % Factor to increase figure size (dashed line hack)
fig_width = 8.6*fac; % in cm
fig_height = 4.5*fac;
fontsize = 8*fac;
lw = .8*fac;

fh=figure(5);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');  
   

 axes('Position',[ml,mb,spanx,spany]);
 hold on
     plot(s1.m,xtop,'-g','linewidth',lw)
     plot(-s1.a1,xtop,'-b','linewidth',lw)
    plot(-s1.A1,xbot,'-r','linewidth',lw)

    ylabel('$x_3$','interpreter','latex','fontsize',fontsize)
%      plot(s1.a1,xtop,':b','linewidth',lw)
%      plot(s1.A1,xbot,':r','linewidth',lw)
    plot(zx,0*zx,'--k','linewidth',lw/2)
    plot(0*zy,zy,'--k','linewidth',lw/2)
    axis([-1 1 -ELL-1 ell+1])
    set(gca,'FontSize',fontsize,...
       'ytick',[-ELL 0 ell],'xtick',[-1 0 1],'ticklabelinterpreter','latex');
       text(-.95,39,['$f=$',num2str(fnew1/10^9,'%.1f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
     text(-.95,29,['$k=$',num2str(kvec(1)/10^9,'%.3f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
%     text(-.95,29,'$(a)$','interpreter','latex','verticalalignment','top',...
%     'horizontalalignment','left','fontsize',10*fac)
   hold off
   legend('$m$','$A_1$','$\overline{A_1}$','interpreter','latex','location','SW')
legend boxoff
    
 axes('Position',[ml+spanx+pr,mb,spanx,spany]);
 hold on
    plot(s2.a3,xtop,'-b','linewidth',lw)
    plot(s2.A3,xbot,'-r','linewidth',lw)
    plot(s2.a1,xtop,':b','linewidth',lw)
    plot(s2.A1,xbot,':r','linewidth',lw)
    plot(s2.m,xtop,'-g','linewidth',lw)
    plot(zx,0*zx,'--k','linewidth',lw/2)
    plot(0*zy,zy,'--k','linewidth',lw/2)
        set(gca,'FontSize',fontsize,...
       'ytick',[],'xtick',[-1 0 1],'ticklabelinterpreter','latex');
    axis([-1 1 -ELL-1 ell+1])
%            text(-.95,-1,['$f=$',num2str(fnew2/10^9,'%.1f')],'interpreter','latex',...
%     'verticalalignment','top','horizontalalignment','left',...
%     'fontsize',7*fac);
%      text(-.95,-9,['$k=$',num2str(kvec(2)/10^9,'%.3f')],'interpreter','latex',...
%     'verticalalignment','top','horizontalalignment','left',...
%     'fontsize',7*fac);
    text(-.95,29,'$(b)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',10*fac)
 hold off
 
  axes('Position',[ml+2*(spanx+pr),mb,spanx,spany]);
 hold on
    plot(s3.a3,xtop,'-b','linewidth',lw)
    plot(s3.A3,xbot,'-r','linewidth',lw)
    plot(s3.a1,xtop,':b','linewidth',lw)
    plot(s3.A1,xbot,':r','linewidth',lw)
    plot(s3.m,xtop,'-g','linewidth',lw)
    plot(zx,0*zx,'--k','linewidth',lw/2)
    plot(0*zy,zy,'--k','linewidth',lw/2)
    set(gca,'FontSize',fontsize,...
       'ytick',[],'xtick',[-1 0 1],'ticklabelinterpreter','latex');
    axis([-1 1 -ELL-1 ell+1])
%     text(.05,-1,['$f=$',num2str(fnew3/10^9,'%.1f')],'interpreter','latex',...
%     'verticalalignment','top','horizontalalignment','left',...
%     'fontsize',7*fac);
%      text(.05,-9,['$k=$',num2str(kvec(3)/10^9,'%.3f')],'interpreter','latex',...
%     'verticalalignment','top','horizontalalignment','left',...
%     'fontsize',7*fac);
    text(-.95,29,'$(c)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',10*fac)

hold off
   
%     text(1,1,txtvec(ii),'interpreter','latex','verticalalignment','top',...
%     'horizontalalignment','left','fontsize',10*fac)
%     text(locvec(ii),1,['$f=$',num2str(fnew/10^9,'%.1f')],'interpreter','latex',...
%     'verticalalignment','top','horizontalalignment','left',...
%     'fontsize',7*fac);
%      text(locvec(ii),.88,['$k=$',num2str(k/10^9,'%.3f')],'interpreter','latex',...
%     'verticalalignment','middle','horizontalalignment','left',...
%     'fontsize',7*fac);


set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
print(fh,'-dpdf','wave_structure_plot_N_12_tom_2.pdf'); 