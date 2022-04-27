%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

N = 6;
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

    kt1 = [.42074 .42609 .43177 .43712 .44147 .44916 .45652 .46321 .46957 .47559];
ft1 = [182.683 187.1504 191.9146 196.3706 199.8996 205.5238 209.816 213.2231 216.3067 219.1702];
 kt2 = [.42074 .42408 .42508 .43043 .43579 .44147 .44783 .45318 .45853 .46388 .46923 .47458];
 ft2 = [194.5097 196.0963 196.5745 199.1497 201.7959 204.7829 208.6479 212.6085 217.0746 221.8059 226.6863 231.6691];
 

mt1 = zeros(size(kt1));
a1 = mt1;
fv1 = mt1;
  
mt2 = zeros(size(kt2));
a2 = mt2;
fv2 = mt2;

parfor ii = 1:length(kt1)
    k = kt1(ii)*10^9;
    f = ft1(ii)*10^9;
    [vout,fnew] = mode_energy_one_layer(f,k,par,N);
    a1(ii) = vout(2)/vout(1);
    fv1(ii)  = fnew/10^9;
end
parfor ii = 1:length(kt2)
    k = kt2(ii)*10^9;
    f = ft2(ii)*10^9;
    [vout,fnew] = mode_energy_one_layer(f,k,par,N);
    a2(ii) = vout(2)/vout(1);
    fv2(ii)  = fnew/10^9;
end
%% Calculate asymptotic solution
knew1 = linspace(min(kt2),max(kt1),200);
anew1 = interp1(kt1,a1,knew1);

knew2 = linspace(min(kt2),max(kt2),200);
anew2 = interp1(kt2,a2,knew2);

K = par.cs/par.beta;
Om = par.cs*K;
H = par.g*par.h0*par.beta/par.cs^2;
G = par.cl/par.cs;
an = 0;
kstar = sqrt(-H+1/2*(1-2*an^2+sqrt(1-4*H+4*an^2*(G^2-1))));
wstar = sqrt(kstar^2+G^2*an^2);

kplot1 = knew1*10^9;
kplot2 = knew2*10^9;

del1 = (kplot1/K-kstar)/par.eps;
f1= -1./(1*(del1-2*del1*wstar-sqrt(del1.^2.*(2*wstar-1).^2+wstar)));

del2 = (kplot2/K-kstar)/par.eps;
f2= 1./(1*(del2-2*del2*wstar+sqrt(del2.^2.*(2*wstar-1).^2+wstar)));



%% Plotting

   %Figure specs
ml = 0.06; % Margin left
mr = 0.00; % Margin right
mu = 0.06; % Margin top
md = 0.17;  % Margin bottom
pb = 0.08; % Interaxes padding bottom
pr = 0.01; % Interaxes padding right
Nplots = 1;
spanx = (1-ml-mr-(Nplots-1)*pr)/Nplots;
spany = (1-mu-md-(1-1)*pb)/1;
fac = 4; % Factor to increase figure size (dashed line hack)
fig_width = 7.6*fac; % in cm
fig_height = 4.5*fac;
fontsize = 8*fac;
lw = 1.5*fac;

fh=figure(2);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
   


hold on
plot(knew1,real(anew1),'-b','linewidth',lw,'handlevisibility','off')
plot(knew2,real(anew2),'-b','linewidth',lw,'handlevisibility','off')
plot(knew1,f1,'-.r','linewidth',lw,'handlevisibility','off')
plot(knew2,f2,'-.r','linewidth',lw,'handlevisibility','off')


plot(nan,nan,'-b','linewidth',lw)
plot(nan,nan,'-r','linewidth',lw)
hold off
 %axis([min([knew1(1) knew2(1)]) max([knew1(end) knew2(end)]) 0 6])
 axis([.420 .475 0 6])

    xlabel('Wavenumber (rad/nm)','interpreter','latex','fontsize',fontsize)
    ylabel('$||A_0||/||m_0||$','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');%,'ytick',[0 .5 1],'xtick',[210 220 230]);
    set(gca,'fontsize',fontsize)



% ov = 0:.05:1;
% a1 = 209.5*ones(size(ov));
% a2 = 227.4*ones(size(ov));
% hold on
% plot(a1,ov,'--k','linewidth',lw/2,'handlevisibility','off')
% plot(a2,ov,'--k','linewidth',lw/2,'handlevisibility','off')
% hold off


% text(227.4,1,'anticrossing','interpreter',...
%       'latex',...
%       'verticalalignment','middle','horizontalalignment','center',...
%       'fontsize',fontsize,'color','black');
%   text(209.5,1,'anticrossing','interpreter',...
%       'latex',...
%       'verticalalignment','middle','horizontalalignment','center',...
%       'fontsize',fontsize,'color','black');
%   text(213,.863,'$m$','interpreter',...
%       'latex',...
%       'verticalalignment','bottom','horizontalalignment','left',...
%       'fontsize',fontsize,'color','black');
% 
%   text(213,.44,'$A_3$','interpreter',...
%       'latex',...
%       'verticalalignment','bottom','horizontalalignment','left',...
%       'fontsize',fontsize,'color','black');
%   
%     text(213,.05,'$\overline{A_3}$','interpreter',...
%       'latex',...
%       'verticalalignment','bottom','horizontalalignment','left',...
%       'fontsize',fontsize,'color','black');
  
legend('SCM','Asymptotics','location','north','interpreter','latex')

legend boxoff 
set(gcf,'Resize','off')

print(fh,'-dpdf','anticrossing_resonant_wave_structure.pdf');
