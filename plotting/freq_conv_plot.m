%Plots convergence figure
 
%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

%f_init = 137.2174;
% f_init = 167.7488;
% f_init = 300.5974;
% f_init = 249.8971;
% f_init = 199.5158;
%f_init = 341.3463;

fvec = [137 167.75  200]*10^9;
kvec = [300 300 300]*10^6;
    f1 = figure;
    f1.Renderer = 'Painters';
    
        
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
fig_height = 4.6*fac;
fontsize = 8*fac;
 
set(f1,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
set(f1,'Resize','off')
    hold on


for ii = 1:length(fvec)

f = fvec(ii);
k = kvec(ii);
    load_data = ['conv_data_v2_',num2str(k/10^6,'%0.f'),'_f_',...
    num2str(f/(10^9),'%0.f')];
    data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\',load_data,'.mat'];
    load(data_dir,'Ntest','k','w','k_init','f_init','freq_vec','freq_conv_err');
    freq_conv_diff = mp('0').*Ntest;
    for hh = 1:length(Ntest)
        freq_conv_diff(hh) = abs(freq_vec(end)-freq_vec(hh))/freq_vec(end);
    end
if ii==1   
plot(Ntest,real(log(freq_conv_diff)),'sg','MarkerSize',6*fac)
elseif ii==2
    plot(Ntest,real(log(freq_conv_diff)),'ob','MarkerSize',5*fac)
else
    plot(Ntest,real(log(freq_conv_diff)),'*r','MarkerSize',6*fac)
end

end
ylabel('$\log_{10}$(error)','interpreter','latex')


xlabel('$N$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex','ytick',[-30,-20,-10],...
    'xtick',[5,10,15,20],'FontSize',fontsize)

     axis([5,20,-30,-4]);

set(f1,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
set(f1,'Resize','off')
hold off

leg1 = legend([num2str(fvec(1)/(10^9),'%.0f'),' GHz'],...
    [num2str(fvec(2)/(10^9),'%.0f'),' GHz'],...
    [num2str(fvec(3)/(10^9),'%.0f'),' GHz'],'interpreter','latex',...
    'location','NE');

set(leg1,'Box','off')

print(f1,'-dpdf','freq_conv_plot.pdf'); 
