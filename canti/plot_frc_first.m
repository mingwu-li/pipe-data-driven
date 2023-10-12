function plot_frc_first(FRCmain,FRCisola,Collmain,Collisola)

thm =struct();
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};
% color
ST = cell(2,1);
ST{1} = {'b--','LineWidth',1.5}; % unstable
ST{2} = {'b-','LineWidth',1.5};  % stable
legs = 'SSM-unstable';
legu = 'SSM-stable'; 

% ssm - main
figure; 
colors = colororder;
% plot(BBCInfo.frequency, BBCInfo.amplitude,'Color',colors(1,:),'Linewidth',2,'DisplayName', 'Backbone - SSMLearn');
ax1 = gca;
FRC = FRCmain;
SNidx = FRC.SNidx;
HBidx = FRC.HBidx;
FRC.st = double(FRC.st);
FRC.st(HBidx) = nan;
FRC.st(SNidx) = nan;
hold(ax1,'on');
plot_stab_lines(FRC.om,FRC.Aout(:,1),FRC.st,ST); %,legs,legu
SNfig = plot(FRC.om(SNidx),FRC.Aout(SNidx,1),thm.SN{:});
set(get(get(SNfig,'Annotation'),'LegendInformation'),...
'IconDisplayStyle','off');
HBfig = plot(FRC.om(HBidx),FRC.Aout(HBidx,1),thm.HB{:});
set(get(get(HBfig,'Annotation'),'LegendInformation'),...
'IconDisplayStyle','off');   
% ssm - isola
FRC = FRCisola;
SNidx = FRC.SNidx;
HBidx = FRC.HBidx;
FRC.st = double(FRC.st);
FRC.st(HBidx) = nan;
FRC.st(SNidx) = nan;
hold(ax1,'on');
plot_stab_lines(FRC.om,FRC.Aout(:,1),FRC.st,ST,legs,legu);
SNfig = plot(FRC.om(SNidx),FRC.Aout(SNidx,1),thm.SN{:});
set(get(get(SNfig,'Annotation'),'LegendInformation'),...
'IconDisplayStyle','off');
HBfig = plot(FRC.om(HBidx),FRC.Aout(HBidx,1),thm.HB{:});
set(get(get(HBfig,'Annotation'),'LegendInformation'),...
'IconDisplayStyle','off');   
xlabel('$\Omega$','Interpreter','latex'); 
ylabel('$||u_1||_{\infty}$','Interpreter','latex'); 
set(gca,'FontSize',14);
grid on; axis tight; 

% coll-main and isola
h1 = plot(Collmain.omega,Collmain.yend,'ms','DisplayName','COCO-unstable');
om = Collisola.omega; 
am = Collisola.yend; st = Collisola.stabs;
plot(om(st),am(st),'ro','DisplayName','COCO-stable');
plot(om(~st),am(~st),'ms','DisplayName','COCO-unstable');
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';


end