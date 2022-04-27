%Generates figure with relative amplitudes of elements of anticrossing

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

N = 12;

%Elastic parameters for the top material
par.cs = mp('3800'); %transverse wave speed m/sec
par.cl = mp('7200'); %Longitudinal wave Speed m/sec
par.rho = mp('5170'); %Density kg/m^3
par.ell = mp('30*10^-9');
par.c11 = par.rho*par.cl^2;
par.c12 = par.rho*(par.cl^2-2*par.cs^2);
par.c44 = par.rho*par.cs^2;

%Magnetic parameters for the top material
par.h0 = mp('.25');     %applied field
par.ms = mp('140000'); %saturation magnetization
par.lam = mp('3')*mp('10^-16'); %exchange length squared
par.mu0 = 4*mp('pi')*mp('10^-7'); %Vacuum permeability
par.g = mp('1.760859644')*mp('10^11'); %gyromagnetic ratio of an electron
par.b2 = mp('7000000'); %Magnetoelastic coupling constant


%Physical parameters for the bottom material
par.CL = mp('6400');
par.CS = mp('3500');
par.RH = mp('7080');
par.C11 = par.RH*par.CL^2;
par.C12 = par.RH*(par.CL^2-2*par.CS^2);
par.C44 = par.RH*par.CS^2;
par.ELL = mp('50*10^-9');

    kt = [.3602 .36161 .36342 .36422 .36583 .36744 .36905 .36945 .36965 .37085 .37106 .37226 .37387 .37548...
        .37709 .37869 .3803 .38191 .38352 .38513 .38673 .38754 .38834...
        .38995 .39176 .39337 .39477 .39538 .39678 .39759 .39859];
ft = [202.5804 203.4323 204.5457 205.0493 206.0791 207.1494 208.2757 208.5665 208.7124 209.659 209.8161 210.8005 212.2058...
        213.6856 215.3376 216.959 218.6342 220.3384 222.0562 223.7639...
        225.3869 226.1012 226.7138 227.7327 228.7608 229.6445 230.4089 230.735 231.4939 231.9265 232.4667];
 
 

mt = zeros(size(kt));
a1t = mt;
a2t = mt;
a3t = mt;
A1t = mt;
A2t = mt;
A3t = mt;

fv = mt;
    
parfor ii = 1:length(kt)
    k = kt(ii)*10^9;
    f = ft(ii)*10^9;
    [vout,fnew] = mode_structure(f,k,par,N);
    mt(ii) = vout(1);
    a1t(ii) = vout(2);
    a2t(ii) = vout(3);
    a3t(ii) = vout(4);
    A1t(ii) = vout(5);
    A2t(ii) = vout(6);
    A3t(ii) = vout(7);
    fv(ii)  = fnew/10^9;
end



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
lw = .8*fac;

fh=figure(2);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
   

fnew = linspace(min(fv),max(fv),200);
mtnew = interp1(fv,mt,fnew);
a3tnew = interp1(fv,a3t,fnew);
A3tnew = interp1(fv,A3t,fnew);


hold on
plot(fnew,real(mtnew),'-g','linewidth',lw,'handlevisibility','off')
plot(fnew,real(a3tnew),'-b','linewidth',lw,'handlevisibility','off')
plot(fnew,real(A3tnew),'-r','linewidth',lw,'handlevisibility','off')

% plot(kt,real(mt),'-g','linewidth',lw,'handlevisibility','off')
% plot(kt,real(a3t),'-b','linewidth',lw,'handlevisibility','off')
% plot(kt,real(A3t),'-r','linewidth',lw,'handlevisibility','off')

plot(nan,nan,'-g','linewidth',2)
plot(nan,nan,'-b','linewidth',2)
plot(nan,nan,'-r','linewidth',2)
hold off
axis([min(fv) max(fv) 0 1])

    xlabel('Frequency [GHz]','interpreter','latex','fontsize',fontsize)
    ylabel('$||\cdot||$','interpreter','latex','fontsize',fontsize)
    set(gca,'TickLabelInterpreter','latex','ytick',[0 .5 1],'xtick',[210 220 230]);
    set(gca,'fontsize',fontsize)



ov = 0:.05:1;
a1 = 209.5*ones(size(ov));
a2 = 227.4*ones(size(ov));
hold on
plot(a1,ov,'--k','linewidth',lw/2,'handlevisibility','off')
plot(a2,ov,'--k','linewidth',lw/2,'handlevisibility','off')
hold off


text(227.4,1,'anticrossing','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','center',...
      'fontsize',fontsize,'color','black');
  text(209.5,1,'anticrossing','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','center',...
      'fontsize',fontsize,'color','black');
  text(213,.863,'$m$','interpreter',...
      'latex',...
      'verticalalignment','bottom','horizontalalignment','left',...
      'fontsize',fontsize,'color','black');

  text(213,.44,'$A_3$','interpreter',...
      'latex',...
      'verticalalignment','bottom','horizontalalignment','left',...
      'fontsize',fontsize,'color','black');
  
    text(213,.05,'$\overline{A_3}$','interpreter',...
      'latex',...
      'verticalalignment','bottom','horizontalalignment','left',...
      'fontsize',fontsize,'color','black');
  
% legend('$m$ ','$A_3$ ','$\overline{A}_3$',...,
%     'location','eastoutside','interpreter','latex')

% legend boxoff 
%set(gcf,'Resize','off')

print(fh,'-dpdf','anticrossing_zoomed_dip.pdf');