%Generate bar graph for given frequencies and wavenumbers

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34
% 
%      kpts = [.36905 .37709 .38955];
%      fpts = [208.2757 215.3376 227.7327];

kvec    = [mp('.36281*10^9') mp('.37709*10^9') mp('.392*10^9')];
fvec = [mp('204.1718*10^9') mp('215.3376*10^9') mp('228.9*10^9')];
N = 12;

%Elastic parameters for the top material (YIG)
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

%Physical parameters for the bottom material (GGG)
par.CL = mp('6400');
par.CS = mp('3500');
par.RH = mp('7080');
par.C11 = par.RH*par.CL^2;
par.C12 = par.RH*(par.CL^2-2*par.CS^2);
par.C44 = par.RH*par.CS^2;
par.ELL = mp('50*10^-9');

    %Figure specs
ml = 0.06; % Margin left
mr = 0.03; % Margin right
mt = 0.04; % Margin top
mb = 0.17;  % Margin bottom
pb = 0.08; % Interaxes padding bottom
pr = 0.01; % Interaxes padding right
Nplots = length(fvec);
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


str = {'$m$','$A_1$','$A_2$','$A_3$','$\overline{A_1}$','$\overline{A_2}$','$\overline{A_3}$'};

locvec = [4,4,4];
txtvec = {'$(a)$','$(b)$','$(c)$'};

for ii = 1:Nplots
    freq = fvec(ii);
    k = kvec(ii);
    [struc,fnew] = mode_structure(freq,k,par,N);
    axes('Position',[ml+(ii-1)*(spanx+pr),mb,spanx,spany]);
   hold on
  b1 =   bar(1,struc(1));
  b2 =   bar(2:4,struc(2:4));
  b3 =   bar(5:7,struc(5:7));
  set(b1,'FaceColor','g')
  set(b2,'FaceColor','b')
  set(b3,'FaceColor','r')
    if ii==1
        ylabel('$||\cdot||$','interpreter','latex','fontsize',6*fac)
    end
set(gca, 'XTickLabel',str, 'XTick',1:numel(str),'ytick',[],...
    'TickLabelInterpreter','Latex','fontsize',8*fac)

axis([.5 7.5 0 1])
    
    text(1,1,txtvec(ii),'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',10*fac)
    text(locvec(ii),1,['$f=$',num2str(fnew/10^9,'%.1f')],'interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
    'fontsize',7*fac);
     text(locvec(ii),.88,['$k=$',num2str(k/10^9,'%.3f')],'interpreter','latex',...
    'verticalalignment','middle','horizontalalignment','left',...
    'fontsize',7*fac);
end

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
print(fh,'-dpdf','mode_energy_plot_N_12_1.pdf'); 