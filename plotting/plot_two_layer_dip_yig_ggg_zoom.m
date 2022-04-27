 %% plot results 

     data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\double_layer_dip_yig_ggg_N_12_zoom.mat'];
    load(data_dir,'par','wtop1','wtop2','wtop3','wbot1','wbot2','wbot3','wmag','kplot','N');

    kt = [.3602 .36161 .36342 .36422 .36583 .36744 .36905 .36945 .36965 .37085 .37106 .37226 .37387 .37548...
        .37709 .37869 .3803 .38191 .38352 .38513 .38673 .38754 .38834...
        .38995 .39176 .39337 .39477 .39538 .39678 .39759 .39859];
ft = [202.5804 203.4323 204.5457 205.0493 206.0791 207.1494 208.2757 208.5665 208.7124 209.659 209.8161 210.8005 212.2058...
        213.6856 215.3376 216.959 218.6342 220.3384 222.0562 223.7639...
        225.3869 226.1012 226.7138 227.7327 228.7608 229.6445 230.4089 230.735 231.4939 231.9265 232.4667];
    
kpts    = [.36281 .37709 .392];
fpts = [204.1718 215.3376 228.9];
 
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

x1 = .35875;
x2 = .36775;
x3 = .36722;
x4 = x1+(x3-x2);
y1 = 216;
y2 = 223;
y3 = 227;
y4 = y1+(y3-y2);

axes('Position',[ml,mb,spanx,spany]);
    
    hold on

    plot(kplot/10^9,sort(real(wbot3))/(2*pi*10^9),'or',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wtop3))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');
     plot(kpts,fpts,'*k','MarkerSize',6*lw,'HandleVisibility','off')
%     rect = plot([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],'-k','linewidth',lw/2,...
%             'handlevisibility','off');
  %  rotate(rect,[0 0 1],10)
 %   rectangle('Position',[.355 215 .012  10],'linewidth',lw/2)
    hold off
    % labels
    axis xy
    xlabel('Wavenumber [rad/nm]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');

 axis([3.61*10^8 3.94*10^8 198*10^9 233*10^9]/10^9)
   set(gca,'ytick',[200,210,220,230,240,250,270,290],...
        'xtick',[0,.36,.37,.38,.39,.4,.41,.42,.43,.44,.45,.6]);  
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'linewidth',lw);     
set(gcf,'Resize','off')




%%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw);
plot(nan,nan,'-r','LineWidth',lw);
plot(nan,nan,'-g','LineWidth',lw);
hold off

hold on
plot(kt,ft,'--k','linewidth',lw,'handlevisibility','off')
hold off

legend('Top Shear','Bottom Shear','Magnetic',...
    'interpreter','latex','location','northoutside','fontsize',fontsize,'orientation','horizontal')
legend boxoff 
set(gcf,'Resize','off')

%Add texts
kpts    = [.36281 .37709 .392];
fpts = [204.1718 215.3376 228.9];
% 
    text(.36281,204.1718 ,'$(a)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
    text(.37709,215.3376,'$(b)$','interpreter','latex','verticalalignment','bottom',...
    'horizontalalignment','right','fontsize',7*fac)
    text(.392,228.9,'$(c)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
%     text(.36355,222.7118,'$e$','interpreter','latex','verticalalignment','bottom',...
%     'horizontalalignment','left','fontsize',7*fac)
%     text(.36873,227.5918,'$f$','interpreter','latex','verticalalignment','top',...
%     'horizontalalignment','left','fontsize',7*fac)
%     text(.36773,228.6617,'$g$','interpreter','latex','verticalalignment','top',...
%     'horizontalalignment','left','fontsize',7*fac)
%     text(x4,y4,'$(a)$','interpreter','latex','verticalalignment','bottom',...
%     'horizontalalignment','left','fontsize',7*fac)
%      text(.402,290,'$(a) \; B_2 \neq 0$','interpreter',...
%       'latex',...
%       'verticalalignment','top','horizontalalignment','left',...
%       'fontsize',fontsize*1.3,'color','black');

print(fh,'-dpdf','two_layer_dip_yig_ggg_N_12_zoom.pdf');