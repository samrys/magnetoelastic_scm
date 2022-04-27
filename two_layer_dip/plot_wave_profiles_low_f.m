%Generate wave profile plots for given frequencies and wavenumbers

%This generates figure for low frequency interactions

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

kvec    = [mp('4.5732*10^6') mp('4.5732*10^6') mp('6.0377*10^6') ;...
            mp('6.2134*10^6') mp('6.8577*10^6') mp('7.6192*10^6')];
fvec = [mp('2.617*10^9') mp('2.7087*10^9') mp('3.5511*10^9'); mp('2.5854*10^9') mp('2.893*10^9') mp('3.5646*10^9') ];
N = 16;

%Elastic parameters for the top material (YIG)
par.cs = mp('3800'); %transverse wave speed m/sec
par.cl = mp('7200'); %Longitudinal wave Speed m/sec
par.rho = mp('5170'); %Density kg/m^3
par.ell = mp('200*10^-9');
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


%Physical parameters for the bottom material (GGG)
par.CL = mp('6400');
par.CS = mp('3500');
par.RH = mp('7080');
par.C11 = par.RH*par.CL^2;
par.C12 = par.RH*(par.CL^2-2*par.CS^2);
par.C44 = par.RH*par.CS^2;
par.ELL = mp('300*10^-9');


xtop = linspace(0,par.ell*10^6,200);
xbot = linspace(-par.ELL*10^6,0,200);

for ii = 1:2
    for jj = 1:3
[s(ii,jj),fnew(ii,jj)] = wave_profiles(fvec(ii,jj),kvec(ii,jj),par,N,xtop,xbot);
    end
end

str = {'$(a)$','$(b)$','$(c)$';'$(d)$','$(e)$','$(f)$'};


zx = -1.05:.05:1.05; %x axis
zy = -.31:.5:.21;      % y axis

%% Plotting
    %Figure specs
ml = 0.1; % Margin left
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
fig_height = 8*fac;
fontsize = 8*fac;
lw = .8*fac;

fh=figure(3);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');  

   
xtop = double(xtop);
xbot = double(xbot);
for ii = 1:2
 for jj =1:3
 axes('Position',[ml+(jj-1)*(pr+spanx),mb+(Ny-ii)*(pb+spany),spanx,spany]);
 hold on
    plot(s(ii,jj).a3,xtop,'-b','linewidth',1.6*lw)
    plot(s(ii,jj).A3,xbot,'-r','linewidth',lw*1.6)
    plot(s(ii,jj).a2,xtop,'-.b','linewidth',lw*1.6)
    plot(s(ii,jj).A2,xbot,'-.r','linewidth',lw*1.6)
    plot(s(ii,jj).m,xtop,'-g','linewidth',lw*1.6)
    set(gca,'FontSize',fontsize,...
       'ytick',[min(xbot) 0 max(xtop)],'xtick',[-1 0 1],'ticklabelinterpreter','latex');


if jj==1
    ylabel('$x_3$','interpreter','latex','fontsize',fontsize)
 
 else
    set(gca,'yticklabel','')
 end
 
 if ii==1
     set(gca,'xticklabel','')
     if jj==3
            legend('$A_3$','$\overline{A_3}$','$A_2$','$\overline{A_2}$','$m$','interpreter','latex','location','SE')
            legend boxoff
     end
 end
 
    plot(zx,0*zx,'--k','linewidth',lw/2,'handlevisibility','off')
    plot(0*zy,zy,'--k','linewidth',lw/2,'handlevisibility','off')
    axis([-1 1 -.305 .205])
   text(-.97,.08,['$f=$',num2str(fnew(ii,jj)/10^9,'%.2f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
     text(-.97,.13,['$k=$',num2str(kvec(ii,jj)/10^6,'%.2f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
    text(-.97,.2,str(ii,jj),'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',10*fac)
   hold off
 end
end

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
%print(fh,'-dpdf','wave_profile_plot_N_16_low_f.pdf'); 