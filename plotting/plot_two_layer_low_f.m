 %% plot results 
 %Generates plots for low frequency anticrossings
 
%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34
 
 
 data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\double_layer_dip_yig_ggg_N_16_low_f_v2.mat'];
 load(data_dir,'par','wtop1','wtop2','wtop3','wbot1','wbot2','wbot3','wmag','kplot','N');

 kvec    = [mp('4.5732*10^6') mp('4.5732*10^6') mp('6.0377*10^6') ; mp('6.2134*10^6') mp('6.8577*10^6') mp('7.6192*10^6')];
fvec = [mp('2.617*10^9') mp('2.7087*10^9') mp('3.5511*10^9'); mp('2.5854*10^9') mp('2.893*10^9') mp('3.5646*10^9') ];
 
kpts = double(kvec)/10^6;
fpts = double(fvec)/10^9;

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
    plot(kplot/10^6,sort(real(wtop2))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
        plot(kplot/10^6,sort(real(wtop3))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^6,sort(real(wbot2))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
        plot(kplot/10^6,sort(real(wbot3))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^6,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kpts,fpts,'*k','MarkerSize',6*lw,'HandleVisibility','off')

    
    text(kpts(1,1),fpts(1,1),'$(a)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
    text(kpts(1,2),fpts(1,2),'$(b)$','interpreter','latex','verticalalignment','bottom',...
    'horizontalalignment','right','fontsize',7*fac)
    text(kpts(1,3),fpts(1,3),'$(c)$','interpreter','latex','verticalalignment','bottom',...
    'horizontalalignment','right','fontsize',7*fac)
    text(kpts(2,1),fpts(2,1),'$(d)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
    text(kpts(2,2),fpts(2,2),'$(e)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
    text(kpts(2,3),fpts(2,3),'$(f)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)

    hold off
    % labels
    axis xy
    xlabel('Wavenumber [rad/$\mu$m$^{-1}$]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
%    title('Dispersion Curves','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');


 axis([4 9 2.1 4])
%     set(gca,'ytick',[230,250,270,290],...
%         'xtick',[0,.36,.37,.38,.39,.4,.41,.42,.43,.44,.45,.6]);  
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex');     
set(gcf,'Resize','off')


%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw);
plot(nan,nan,'-r','LineWidth',lw);
plot(nan,nan,'-g','LineWidth',lw);
%plot(nan,nan,'--k','LineWidth',lw);
hold off

legend('Top Shear','Bottom Shear','Magnetic',...
    'interpreter','latex','location','SE','fontsize',fontsize)
legend boxoff 
 set(gcf,'Resize','on')


print(fh,'-dpdf','two_layer_dip_yig_ggg_N_16_low_f');