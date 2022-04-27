
     data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\one_layer_dip_yig_N_16_b2_nonzero_v2.mat'];
     load(data_dir,'par','wlong','wmag','wtran','kplot','N');
 
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

kvec    = [mp('.35304*10^9') mp('.35906*10^9') mp('.37351*10^9'); mp('.36057*10^9') mp('.36147*10^9') mp('.33769*10^9')];
fvec = [mp('196.9032*10^9') mp('194.0625*10^9') mp('211.9213*10^9'); mp('212.7947*10^9') mp('205.0177*10^9') mp('188.1002*10^9')];

kpts = double(kvec/10^9);
fpts = double(fvec/10^9);

fh=figure(2);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

set(gcf,'Resize','off')

axes('Position',[ml,mb,spanx,spany]);
    
   
    hold on


    plot(kplot/10^9,sort(real(wtran))/(2*pi*10^9),'ob',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kplot/10^9,sort(real(wmag))/(2*pi*10^9),'og',...
        'MarkerSize',lw,'HandleVisibility','off');
    plot(kpts,fpts,'*k','MarkerSize',6*lw,'HandleVisibility','off')
    hold off
    % labels
    axis xy
    xlabel('Wavenumber [rad/nm]','interpreter','latex','fontsize',fontsize)
    ylabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex');

   data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data\single_layer_b2_nonzero.mat'];
     save(data_dir,'par','wlong','wtran','wmag','kplot','N');
 axis([3.3*10^8 3.85*10^8 178*10^9 220*10^9]/10^9)    
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'linewidth',lw);     
 set(gca,'xtick',[.32 .34 .36 .38],'ytick',[180 200 220 240])
set(gcf,'Resize','off')

%Create legend
hold on
plot(nan,nan,'-b','LineWidth',lw/2);
plot(nan,nan,'-g','LineWidth',lw/2);
hold off

legend('Shear','Magnetic',...
    'interpreter','latex','location','SE','fontsize',fontsize)
legend boxoff 
set(gcf,'Resize','off')

    text(kpts(1,1),fpts(1,1),'$(a)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
    text(kpts(1,2),fpts(1,2),'$(b)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
    text(kpts(1,3),fpts(1,3),'$(c)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)
    text(kpts(2,1),fpts(2,1),'$(d)$','interpreter','latex','verticalalignment','bottom',...
    'horizontalalignment','right','fontsize',7*fac)
    text(kpts(2,2),fpts(2,2),'$(e)$','interpreter','latex','verticalalignment','bottom',...
    'horizontalalignment','right','fontsize',7*fac)
    text(kpts(2,3),fpts(2,3),'$(f)$','interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','fontsize',7*fac)

print(fh,'-dpdf','one_layer_dip_yig_30_b2_nonzero_n_16_zoom.pdf');
%set(gcf,'Resize','on')
