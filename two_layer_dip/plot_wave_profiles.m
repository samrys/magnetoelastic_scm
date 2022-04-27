%Generate plots for wave profiles for given frequencies and wavenumbers
%Two layers YIG and GGG

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

%wavenumbers (rad/m) and frequencies (Hz) at which to generate the wave profiles
kvec    = [mp('.36281*10^9') mp('.37709*10^9') mp('.392*10^9')];
fvec = [mp('204.1718*10^9') mp('215.3376*10^9') mp('228.9*10^9')];

N = 16;

%Elastic parameters for the top material
par.cs = mp('3800'); %transverse wave speed m/sec
par.cl = mp('7200'); %Longitudinal wave Speed m/sec
par.rho = mp('5170'); %Density kg/m^3
par.ell = mp('30*10^-9'); %Thickness
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
par.ELL = mp('50*10^-9'); %Thickness

xtop = linspace(0,30,150); %Independent var for top material
xbot = linspace(-50,0,300); %Independent var for bottom material

clear s
clear fnew

for ii = 1:3
[s(ii),fnew(ii)] = wave_profiles(fvec(ii),kvec(ii),par,N,xtop,xbot);
end

zx = -1.05:.05:1.05; %x axis
zy = -55:.5:35;      % y axis

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
fig_height = 4.2025*fac;
fontsize = 8*fac;
lw = 1.5*fac;

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');  
   
str = {'$(a)$','$(b)$','$(c)$'};
   
for ii = 1:3
 axes('Position',[ml+(ii-1)*(spanx+pr),mb,spanx,spany]);
 hold on
    plot(s(ii).a3,xtop,'-b','linewidth',lw)
    plot(s(ii).A3,xbot,'-r','linewidth',lw)
    plot(s(ii).m,xtop,'-g','linewidth',lw)
    ylabel('$x_3$','interpreter','latex','fontsize',fontsize)
    plot(zx,0*zx,'--k','linewidth',lw/2)
    plot(0*zy,zy,'--k','linewidth',lw/2)
    axis([-1 1 -51 31])
    set(gca,'FontSize',fontsize,...
       'ytick',[-50 0 30],'xtick',[-1 0 1],'ticklabelinterpreter','latex');
%        text(-.95,-1,['$f=$',num2str(fnew(ii)/10^9,'%.1f')],'interpreter','latex',...
%     'verticalalignment','top','horizontalalignment','left',...
%     'fontsize',7*fac);
%      text(-.95,-9,['$k=$',num2str(kvec(ii)/10^9,'%.3f')],'interpreter','latex',...
%     'verticalalignment','top','horizontalalignment','left',...
%     'fontsize',7*fac);
    text(-.95,29,str(ii),'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',10*fac)
   hold off
   
   if ii ==1
        legend('$A_3$','$\overline{A_3}$','$m$','interpreter','latex','location','SW')
        legend boxoff
   else
       set(gca,'ytick','')
   
   end
end   

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
%print(fh,'-dpdf','wave_profiles_two_layer_yig_ggg_N_16.pdf'); 