function [frcs,color,thm,ST] = plot_isola_and_main(BBCInfo,FRC1st,FRC2nd,idxplot)

thm =struct();
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};
color = {'r','k','m','b','g'};
% color
ST = cell(2,1);
ST{1} = {'b--','LineWidth',1.5}; % unstable
ST{2} = {'b-','LineWidth',1.5};  % stable
legs = 'SSM-unstable';
legu = 'SSM-stable';
figure; 
colors = colororder;
ax1 = gca;
for k=idxplot
    FRC = FRC1st{k};
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
    xlabel('$\Omega$','Interpreter','latex'); 
    ylabel('$||u_1||_{\infty}$','Interpreter','latex'); 
    set(gca,'FontSize',14);
    grid on; axis tight; 
end
%%
% figure; ax1 = gca;
for k=idxplot
    FRC = FRC2nd{k};
    SNidx = FRC.SNidx;
    HBidx = FRC.HBidx;
    FRC.st = double(FRC.st);
    FRC.st(HBidx) = nan;
    FRC.st(SNidx) = nan;
    hold(ax1,'on');
    if k==1
        plot_stab_lines(FRC.om,FRC.Aout(:,1),FRC.st,ST,legs,legu);
    else
        plot_stab_lines(FRC.om,FRC.Aout(:,1),FRC.st,ST);
    end
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
end
plot(BBCInfo.frequency, BBCInfo.amplitude,'Color',colors(1,:),'Linewidth',2,'DisplayName', 'Backbone - SSMLearn');
frcs = gcf;
end


