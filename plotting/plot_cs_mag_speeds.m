    %Figure specs
ml = 0.12; % Margin left
mr = 0.09; % Margin right
mt = 0.1; % Margin top
mb = 0.17;  % Margin bottom
pb = 0.08; % Interaxes padding bottom
pr = 0.01; % Interaxes padding right
Nplots = 1;
spanx = (1-ml-mr-(Nplots-1)*pr)/Nplots;
spany = (1-mt-mb-(1-1)*pb)/1;
fac = 4; % Factor to increase figure size (dashed line hack)
fig_width = 8.6*fac; % in cm
fig_height = 4.5*fac;
fontsize = 8*fac;
lw = 1.1*fac;

clear h;

t = 10^9;

kvec = linspace(0,1*10^9,250);

% Speeds of sound in meters/second
ys            =   3800        ;
ns            =   3000        ;
gs           =   3500       ;
ss           =   4500        ;

%Magnetism
g = 1.7608596e11;
mu0 = 4*pi*10^-7;
l1 = 7.72e-9;
M1 = 4.8e5;
h1 = .65;
om1 = sqrt((g*h1+g.*mu0.*l1^2*M1.*kvec.^2).*(g*h1+g.*mu0.*l1^2*M1.*kvec.^2+g*mu0*M1))/(2*pi);

l2 = 1.73e-08;
M2 = 1.4e5;
h2 = .25;
om2 = sqrt((g*h2+g.*mu0.*l2^2*M2.*kvec.^2).*(g*h2+g.*mu0.*l2^2*M2.*kvec.^2+g*mu0*M2))/(2*pi);


fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');


axes('Position',[ml,mb,spanx,spany]);
wnl = ys*kvec/(2*pi);
wns = ns*kvec/(2*pi);
snl = gs*kvec/(2*pi);
sns = ss*kvec/(2*pi);

hold on

plot(kvec/t,snl/t,'-k','linewidth',lw,'handlevisibility','off');
plot(kvec/t,sns/t,'-k','linewidth',lw,'handlevisibility','off');
plot(kvec/t,wnl/t,'-k','linewidth',lw,'handlevisibility','off');
plot(kvec/t,wns/t,'-k','linewidth',lw,'handlevisibility','off');
plot(kvec/t,om1/t,'-.k','linewidth',lw,'handlevisibility','off');
plot(kvec/t,om2/t,'-.k','linewidth',lw,'handlevisibility','off');

plot(nan,nan,'-k','linewidth',.5*lw)
plot(nan,nan,'-.k','linewidth',.5*lw)

%plot(ones(size(kvec))*1/l1*10^-9,snl/t,'--k','linewidth',lw)
%plot(ones(size(kvec))*1/l2*10^-9,snl/t,'--k','linewidth',lw)

h = legend('Elastic','Magnetic','interpreter','latex','location','SE');
xlabel('Wavenumber (rad/nm)','interpreter','latex')
ylabel('Frequency (GHz)','interpreter','latex')
 axis([.26 .58 100 380])

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

legend box off

text(.58,276,'Ni','interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*fac);
 text(.58,322,'GGG','interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*fac);
 text(.58,350,'YIG','interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*fac);
 text(.58,363,'Ni','interpreter','latex',...
     'verticalalignment','bottom','horizontalalignment','left',...
     'fontsize',7*fac);
 text(.53,380,'Si$_2$N$_3$','interpreter','latex',...
     'verticalalignment','bottom','horizontalalignment','left',...
     'fontsize',7*fac);
  text(.5,380,'YIG','interpreter','latex',...
     'verticalalignment','bottom','horizontalalignment','right',...
     'fontsize',7*fac);
 
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
       'linewidth',fac/2,'ytick',[0,100,200,300,400,500,600],'xtick',[0,.2,.3,.4,.5,.6]);
% 
%    r1 = get(h,'Position');
%    rect = [.12 .1 0 0];
%    set(h,'Position',rect+r1);


% text(500,sl/1000,[num2str(sl),' m/s'],'interpreter','latex',...
%      'verticalalignment','middle','horizontalalignment','left',...
%      'fontsize',10*f);
%  text(500,ss/1000,[num2str(ss),' m/s'],'interpreter','latex',...
%      'verticalalignment','middle','horizontalalignment','left',...
%      'fontsize',10*f);
%  text(500,nl/1000,[num2str(nl),' m/s'],'interpreter','latex',...
%      'verticalalignment','middle','horizontalalignment','left',...
%      'fontsize',10*f);
%  text(500,ns/1000,[num2str(ns),' m/s'],'interpreter','latex',...
%      'verticalalignment','middle','horizontalalignment','left',...
%      'fontsize',10*f); 
print(fh,'-dpdf','speeds_of_sound_mag.pdf');
