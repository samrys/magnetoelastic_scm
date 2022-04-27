%Plots nonresonant intersection with asymptotic prediction

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34 


   %Figure specs
ml = 0.11; % Margin left
mr = 0.03; % Margin right
mt = 0.03; % Margin top
mb = 0.17;  % Margin bottom
pb = 0.08; % Interaxes padding bottom
pr = 0.08; % Interaxes padding right
Nplots = 2;
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


%% Left figure
axes('Position',[ml,mb,spanx,spany]);

data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\single_layer_simp_j_0_n_1.mat'];
load(data_dir,'par','wtran','wmag','kplot','N');

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
    set(gca,'TickLabelInterpreter','latex');%,'xtick',[.42 .44 .46 .48 .5]);

axis([.482*10^9 .49*10^9 239*10^9 245.9*10^9]/10^9)
    set(gca,'ytick',[240,242,244],...
        'xtick',[.485 .49]);
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'linewidth',fac/2);     

K = par.cs/par.beta;
Om = par.cs*K;
H = par.g*par.h0*par.beta/par.cs^2;
G = par.cl/par.cs;
eps = par.eps;

%Plot asymptotic prediction

aj = 0;
an = par.beta*pi/(par.cs*par.ell);
km = kplot/K;
wj0 = H+km.^2+aj.^2;
wn0 = sqrt(km.^2+G^2.*an.^2);
Bn = km.*wn0./(wn0.^2-(H+km.^2+an.^2).^2);
Cj = km./(wj0.*(wj0.^2-km.^2-G^2*aj^2));

wn2 = km.*Bn./(2.*wn0).*(H+km.^2+an.^2);
wj2 = km.*Cj./2;

w1 = wn0+eps^2.*wn2;
w2 = wj0 + eps.^2.*wj2;

fasymp1 = Om.*w1/(2*pi);
fasymp2 = Om.*w2/(2*pi);

%Include bulk dispersion curves
fv1 = Om/(2*pi)*(H+kplot.^2./K^2+aj.^2);
fv2 = Om/(2*pi)*sqrt(kplot.^2./K^2+G^2*an.^2);

hold on
plot(kplot/10^9,fv1/10^9,'k-','Linewidth',.5*lw,'HandleVisibility','off')
plot(kplot/10^9,fv2/10^9,'k-','Linewidth',.5*lw,'HandleVisibility','off')
plot(kplot/10^9,fasymp1/10^9,'r-.','Linewidth',1.5*lw,'HandleVisibility','off')
plot(kplot/10^9,fasymp2/10^9,'r-.','Linewidth',1.5*lw,'HandleVisibility','off')

hold off

    text(.482,245.9,'$(n,j)=(1,0)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',8*fac)



%% Right figure
axes('Position',[ml+spanx+pr,mb,spanx,spany]);

data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\single_layer_simp_j_1_n_0.mat'];
load(data_dir,'par','wtran','wmag','kplot','N');

hold on
plot(kplot/10^9,sort(real(wtran))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');
hold off
    % labels
    axis xy
    xlabel('Wavenumber [rad/nm]','interpreter','latex','fontsize',fontsize)
   set(gca,'TickLabelInterpreter','latex');%,'xtick',[.42 .44 .46 .48 .5]);

axis([.418*10^9 .427*10^9 192.1*10^9 198.5*10^9]/10^9)
    set(gca,'ytick',[193,195,197],...
        'xtick',[.42 .425]);

 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex')
  
%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw);
plot(nan,nan,'-g','LineWidth',lw);
plot(nan,nan,'r-.','Linewidth',lw);
plot(nan,nan,'k-','linewidth',.5*lw);


legend('Shear','Magnetic','Asymptotic','Uncoupled',...
    'interpreter','latex','location','SE','fontsize',fontsize)
legend boxoff 

%Plot asymptotic prediction

aj = 1*par.beta*pi/(par.cs*par.ell);
an = 0*par.beta*pi/(par.cs*par.ell);
km = kplot/K;
wj0 = H+km.^2+aj.^2;
wn0 = sqrt(km.^2+G^2.*an.^2);
Bn = km.*wn0./(wn0.^2-(H+km.^2+an.^2).^2);
Cj = km./(wj0.*(wj0.^2-km.^2-G^2*aj^2));

wn2 = km.*Bn./(2.*wn0).*(H+km.^2+an.^2);
wj2 = km.*Cj./2;

w1 = wn0+eps^2.*wn2;
w2 = wj0 + eps.^2.*wj2;

fasymp1 = Om.*w1/(2*pi);
fasymp2 = Om.*w2/(2*pi);

%Include bulk dispersion curves
fv1 = Om/(2*pi)*(H+kplot.^2./K^2+aj.^2);
fv2 = Om/(2*pi)*sqrt(kplot.^2./K^2+G^2*an.^2);

hold on
plot(kplot/10^9,fv1/10^9,'k-','Linewidth',.5*lw,'HandleVisibility','off')
plot(kplot/10^9,fv2/10^9,'k-','Linewidth',.5*lw,'HandleVisibility','off')
plot(kplot/10^9,fasymp1/10^9,'r-.','Linewidth',1.5*lw,'HandleVisibility','off')
plot(kplot/10^9,fasymp2/10^9,'r-.','Linewidth',1.5*lw,'HandleVisibility','off');

hold off


% 
    text(.418,198.5,'$(n,j)=(0,1)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',8*fac)

%print(fh,'-dpdf','nonresonant_crossings.pdf');

%set(gcf,'Resize','on')

