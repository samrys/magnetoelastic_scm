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
t = 10^9;

kvec = linspace(0,1*10^9,300);

% Speeds of sound in meters/second
ns            =   2900        ;

%Magnetism
g = 1.7608596e11;
mu0 = 4*pi*10^-7;
l = 7.72e-9;
M = 4.8e5;
b2 = 10000000*2;
rh = 8900;
h0 = .65;
w0 = h0*g;
wM = g*mu0*M;

wm = sqrt((w0+wM*kvec.^2*l.^2).*(w0+wM*(kvec.^2*l.^2-1)));
wl = ns*kvec;
wc = g.*b2.^2.*kvec.^2./(rh.*M).*((w0+wM*kvec.^2*l.^2));
wms = wm.^2;
wls = wl.^2;


w1s = .5.*(wls+wms+sqrt(4*wc+wls.^2-2.*wls.*wms+wms.^2));
w2s = .5.*(wls+wms-sqrt(4*wc+wls.^2-2.*wls.*wms+wms.^2));

w1 = sqrt(w1s);
w2 = sqrt(w2s);


fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');


axes('Position',[ml,mb,spanx,spany]);


vv = linspace(0,max(wl),100);


hold on

plot(kvec/t,w1/(2*pi*t),'-r','linewidth',lw,'handlevisibility','off');
plot(kvec/t,w2/(2*pi*t),'-r','linewidth',lw,'handlevisibility','off');
plot(kvec/t,wl/(2*pi*t),'-.k','linewidth',.5*lw,'handlevisibility','off');
plot(kvec/t,wm/(2*pi*t),'-.k','linewidth',.5*lw,'handlevisibility','off');

plot(nan,nan,'-.k','linewidth',.5*lw);
plot(nan,nan,'-r','linewidth',.5*lw);


%plot(ones(size(kvec))*1/l1*10^-9,snl/t,'--k','linewidth',lw)
%plot(ones(size(kvec))*1/l2*10^-9,snl/t,'--k','linewidth',lw)

h = legend('Uncoupled','Coupled','interpreter','latex','location','SE');
xlabel('Wavenumber [rad/nm]','interpreter','latex')
ylabel('Frequency [GHz]','interpreter','latex')
 axis([0 .5 0 240])
legend boxoff 
text(.65,285,'Ni','interpreter','latex',...
     'verticalalignment','bottom','horizontalalignment','left',...
     'fontsize',7*fac);

 text(.62,400,'Ni','interpreter','latex',...
     'verticalalignment','bottom','horizontalalignment','left',...
     'fontsize',7*fac);

 
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
       'linewidth',fac/2,'ytick',[0,50,100,150,200,250,300,400,500,600],'xtick',[0,.1,.2,.3,.4,.5,.6]);
kl = 0;
 kh = .04;
 wb = 0;
 wh = 13;
 rectangle('Position',[kl wb kh-kl wh-wb],'linewidth',3)
 
axes('Position',[.2,.6,.28,.25]);
box on
hold on
plot(kvec/t,w1/(2*pi*t),'-r','linewidth',lw,'handlevisibility','off');
plot(kvec/t,w2/(2*pi*t),'-r','linewidth',lw,'handlevisibility','off');
plot(kvec/t,wl/(2*pi*t),'-.k','linewidth',.5*lw,'handlevisibility','off');
plot(kvec/t,wm/(2*pi*t),'-.k','linewidth',.5*lw,'handlevisibility','off');
 set(gca,'FontSize',.8*fontsize,'TickLabelInterpreter','latex',...
       'linewidth',fac/2,'ytick',[0,10],'xtick',[0,.02,.04,.05,.6]);
hold off
axis([kl kh wb wh])

% annotation('line',[.12 .2],[.21 .6],'linewidth',fac/2)
% annotation('line',[.18 .45],[.21 .6],'linewidth',fac/2)

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');




print(fh,'-dpdf','nickel_regimes.pdf');
