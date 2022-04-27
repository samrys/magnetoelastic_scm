 %% plot results 

      data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\double_layer_dip_ni_sn_N_12.mat'];
     load(data_dir,'par','wlong_top','wlong_bot','wshear_top','wshear_bot','wmag','kplot','N'); 
 
 
    %Figure specs
ml = 0.11; % Margin left
mr = 0.03; % Margin right
mt = 0.00; % Margin top
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

    plot(kplot/10^9,sort(real(wshear_bot))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wshear_top))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');

    hold off
    % labels
    axis xy
    xlabel('Wavenumber [rad/nm]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');

 axis([3.35*10^8 3.97*10^8 184*10^9 238*10^9]/10^9)
%    set(gca,'ytick',[200,220,240,260,280],...
%         'xtick',[0,.36,.37,.38,.39,.4,.41,.42,.43,.44,.45,.6]);  
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'linewidth',lw);     
set(gcf,'Resize','off')



%%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw);
plot(nan,nan,'-r','LineWidth',lw);
plot(nan,nan,'-g','LineWidth',lw);
hold off

legend('Top Shear','Bottom Shear','Magnetic',...
    'interpreter','latex','location','northoutside',...
    'fontsize',fontsize,'orientation','horizontal')
legend boxoff 
set(gcf,'Resize','off')

%      text(.402,290,'$(a) \; B_2 \neq 0$','interpreter',...
%       'latex',...
%       'verticalalignment','top','horizontalalignment','left',...
%       'fontsize',fontsize*1.3,'color','black');

print(fh,'-dpdf','two_layer_dip_ni_sn_N_12.pdf');