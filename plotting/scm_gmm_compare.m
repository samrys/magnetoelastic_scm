  


%% plot results map

 load('C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\gmm_disp_curves.mat')
    % Need to have ElasticMatrix code on the path for this to work  

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
f = 4; % Factor to increase figure size (dashed line hack)
fig_width = 8.6*f; % in cm
fig_height = 4.5*f;
fontsize = 8*f;

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

set(gcf,'Resize','off')

axes('Position',[ml,mb,spanx,spany]);
    
         hold on
     for ii = 1:58
     f1 = my_model.dispersion_curves(ii).f/10^9;
     k1 = my_model.dispersion_curves(ii).k/10^9;
     plot(k1,f1,'-k','handlevisibility','off','linewidth',.7*f);
     end
     
     plot(nan,nan,'-k')
     plot(nan,nan,'--r')
    
 load('C:\Users\samry\Documents\MATLAB\Elastic\data\no_mag_N_24_compare.mat','kplot','wplot');
   
  %  kplot = double(kv);
    plot(kplot/10^9,sort(wplot)/10^9/(2*pi),'--r','handlevisibility','off','linewidth',f*.7);
    % labels
    axis xy
    xlabel('Wavenumber [rad/nm]','interpreter','latex')
    ylabel('Frequency [GHz]','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex');

%    data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\data\stor_',num2str(n),'.mat'];
%    save(data_dir,'par','kv','omv','A','n');
     axis([1.5*10^8 2.5*10^8 200*10^9 300*10^9]/10^9)
      set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
       'linewidth',f/2,'ytick',[0,100,200,300,400,500,600],...
       'xtick',[0,.1,.15,.2,.25,.3,.4,.5,.6]); 
     
      %% Compare to bulk speeds
nl            =   5500        ;
ns            =   2900        ;

sl           =   7700       ;
ss           =   4600        ;

wnl = nl*kplot/(2*pi);
wns = ns*kplot/(2*pi);
snl = sl*kplot/(2*pi);
sns = ss*kplot/(2*pi);

hold on
%  plot(kplot/10^9,wnl/10^9,'-.g','linewidth',.5*f);
%  plot(kplot/10^9,wns/10^9,'-.g');
  plot(kplot/10^9,snl/10^9,'-.g','linewidth',.7*f);
%  plot(kplot/10^9,sns/10^9,'-.g');
hold off
     
     %%

    
    legend('GMM','SCM','Si$_3$N$_4$ $c_L$','interpreter','latex','location','SE')
    
    print(fh,'-dpdf','scm_gmm_compare.pdf');
