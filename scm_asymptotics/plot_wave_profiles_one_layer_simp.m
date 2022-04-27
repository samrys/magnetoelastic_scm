%Generate plots for given frequencies and wavenumbers

%Compares nonresonant asymptotics with nonresonant wave structures

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

kvec    = [mp('.42241*10^9') mp('.42241*10^9')];
fvec = [mp('195.3017*10^9') mp('195.2672*10^9')];
N = 8;

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

xtop = linspace(0,30,150);

for jj = 1:2
[s(jj),fnew(jj)] = wave_structure_one_layer_simp(fvec(jj),kvec(jj),par,N,xtop);
end



str = {'$(a)$ $\,n=0$','$(b)$ $\,j=1$','$(c)$'};

% locvec = [4,4,4];
% txtvec = {'$(a)$','$(b)$','$(c)$'};

zx = -1.05:.05:1.05;
zy = -5:.5:35;
K = par.cs/par.beta;
Om = par.cs*K;
H = par.g*par.h0*par.beta/par.cs^2;
G = par.cl/par.cs;
eps = par.eps;


%Plot asymptotic prediction

aj = 1*par.beta*pi/(par.cs*par.ell);
an = 0*par.beta*pi/(par.cs*par.ell);
km = kvec(1)/K; %Okay since k values are the same;
wj0 = H+km.^2+aj.^2;
wn0 = sqrt(km.^2+G^2.*an.^2);

%Magnetic correction to elastic bulk wave
Cj = km./(wj0.*(wj0.^2-km.^2-G^2*aj^2));

%Elastic correction to magnetic bulk wave
Bn = km.*wn0./(wn0.^2-(H+km.^2+an.^2).^2);

%Nondimensionalize the vertical variable
znd = xtop*pi/30;

mas = zeros(2,length(znd));
Aas = mas;

mas(1,:) = Bn*eps.*ones(size(znd));
Aas(1,:) = ones(size(znd));
mas(2,:) = -cos(znd);
Aas(2,:) = Cj*eps*cos(znd);

%% Plotting
    %Figure specs
ml = 0.09; % Margin left
mr = 0.02; % Margin right
mt = 0.02; % Margin top
mb = 0.1;  % Margin bottom
pb = 0.05; % Interaxes padding bottom
pr = 0.03; % Interaxes padding right
Nx = 2;
Ny = 1;
spanx = (1-ml-mr-(Nx-1)*pr)/Nx;
spany = (1-mt-mb-(Ny-1)*pb)/Ny;
fac = 4; % Factor to increase figure size (dashed line hack)
fig_width = 8.6*fac; % in cm
fig_height = 4.5*fac;
fontsize = 8*fac;
lw = 1.5*fac;

fh=figure(3);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');  
   
for jj = 1:2
 axes('Position',[ml+(jj-1)*(pr+spanx),mb,spanx,spany]);
 hold on
    plot(s(jj).a3,xtop,'-b','linewidth',lw,'handlevisibility','off')
    plot(s(jj).m2,xtop,'-g','linewidth',lw,'handlevisibility','off')
    plot(mas(jj,:),xtop,'-.r','linewidth',lw,'handlevisibility','off');
    plot(Aas(jj,:),xtop,'-.r','linewidth',lw,'handlevisibility','off');
    set(gca,'FontSize',fontsize,'ytick',[0 30],'xtick',[-1 0 1],'ticklabelinterpreter','latex');

if jj==1
    ylabel('$z$','interpreter','latex','fontsize',fontsize)
 else
    set(gca,'yticklabel','')
 end
 
% plot(zx,0*zx,'--k','linewidth',lw/2)
     plot(0*zy,zy,'--k','linewidth',lw/2,'handlevisibility','off')
    axis([-1.05 1.05 0 31])
  
   text(-1,26,['$f=$',num2str(fnew(jj)/10^9,'%.1f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
     text(-1,23,['$k=$',num2str(kvec(jj)/10^9,'%.3f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
    text(-1,31,str(jj),'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',10*fac)
   hold off
end
hold on
plot(nan,nan,'-b','linewidth',lw)
plot(nan,nan,'-g','linewidth',lw)
plot(nan,nan,'-r','linewidth',lw)

legend('Shear','Magnetic','Asymptotics','location','SE','interpreter','latex')
legend boxoff
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
print(fh,'-dpdf','wave_profile_plot_N_8_asymptotic_comparison.pdf'); 