%Generate plots for given frequencies and wavenumbers for one layer
%material

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

%Vectors at which to determine wave profiles (in rad/nm and Hz)
kvec    = [mp('.35304*10^9') mp('.35906*10^9') mp('.37351*10^9');...
            mp('.36057*10^9') mp('.36147*10^9') mp('.33769*10^9')];
        
fvec = [mp('196.9032*10^9') mp('194.0625*10^9') mp('211.9213*10^9');...
            mp('212.7947*10^9') mp('205.0177*10^9') mp('188.1002*10^9')];
N = 16;

%Elastic parameters for the top material (currently YIG)
par.cs = mp('3800'); %transverse wave speed m/sec
par.cl = mp('7200'); %Longitudinal wave Speed m/sec
par.rho = mp('5170'); %Density kg/m^3
par.ell = mp('30*10^-9'); %Material thickness (paper calls this d)
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

xtop = linspace(0,30,150);

for ii = 1:2
    for jj = 1:3
[s(ii,jj),fnew(ii,jj)] = wave_profiles_one_layer(fvec(ii,jj),kvec(ii,jj),par,N,xtop);
    end
end


str = {'$(a)$','$(b)$','$(c)$';'$(d)$','$(e)$','$(f)$'};

zx = -1.05:.05:1.05;
zy = -5:.5:35;

%% Plotting
    %Figure specs
ml = 0.09; % Margin left
mr = 0.02; % Margin right
mt = 0.02; % Margin top
mb = 0.1;  % Margin bottom
pb = 0.05; % Interaxes padding bottom
pr = 0.03; % Interaxes padding right
Nx = 3;
Ny = 2;
spanx = (1-ml-mr-(Nx-1)*pr)/Nx;
spany = (1-mt-mb-(Ny-1)*pb)/Ny;
fac = 4; % Factor to increase figure size (dashed line hack)
fig_width = 8.6*fac; % in cm
fig_height = 4.5*fac;
fontsize = 8*fac;
lw = .8*fac;

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');  
   
for ii = 1:2
    for jj = 1:3
 axes('Position',[ml+(jj-1)*(pr+spanx),mb+(Ny-ii)*(pb+spany),spanx,spany]);
 hold on
    plot(s(ii,jj).a3,xtop,'-b','linewidth',lw)
    plot(s(ii,jj).m,xtop,'-g','linewidth',lw)
    set(gca,'FontSize',fontsize,'ytick',[0 30],'xtick',[-1 0 1],'ticklabelinterpreter','latex');

if jj==1
    ylabel('$x_3$','interpreter','latex','fontsize',fontsize)
 else
    set(gca,'yticklabel','')
 end
 
 if ii==1
     set(gca,'xticklabel','')
 end
 
     plot(0*zy,zy,'--k','linewidth',lw/2)
    axis([-1.05 1.05 0 31])
  
   text(-1,24,['$f=$',num2str(fnew(ii,jj)/10^9,'%.1f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
     text(-1,19,['$k=$',num2str(kvec(ii,jj)/10^9,'%.3f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
    text(-1,31,str(ii,jj),'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',10*fac)
   hold off
    end
end

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
%print(fh,'-dpdf','wave_profiles_plot_N_16_one_layer.pdf'); 