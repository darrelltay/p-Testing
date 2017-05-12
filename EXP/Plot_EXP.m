function Plot_EXP(Result,XLabelName)
%%%% Result is the result obtained from PTesting_PL
LabelFontsize = 20;
Fontsize = 20;
TextFontsize = 15;
YLabel = 'CDF';

%% Extract Data
Sample = Result.Sample;
%%%% Unadjusted
Xmin_Unadj = Result.Unadjusted.Parameter.Xmin;
Beta_Unadj = Result.Unadjusted.Parameter.Beta;
KS_CSN_Unadj = Result.Unadjusted.PValue.CSN_KS;
KS_IF_Unadj = Result.Unadjusted.PValue.IF_KS;
DN_CSN_Unadj = Result.Unadjusted.PValue.CSN_DN;
DN_IF_Unadj = Result.Unadjusted.PValue.IF_DN;
%%%% Adjusted 
Xmin_Adj = Result.Adjusted.Parameter.Xmin;
Beta_Adj = Result.Adjusted.Parameter.Beta;
KS_CSN_Adj = Result.Adjusted.PValue.CSN_KS;
KS_IF_Adj = Result.Adjusted.PValue.IF_KS;
DN_CSN_Adj = Result.Adjusted.PValue.CSN_DN;
DN_IF_Adj = Result.Adjusted.PValue.IF_DN;
%% Prepare for Plot
Sample = sort(Sample);
EmpCDF = linspace(1,0,size(Sample,1));
%%%% Unadjusted
WoSS = linspace(Xmin_Unadj,max(Sample),50);
WoCDF = EmpCDF(sum(Sample<Xmin_Unadj))*exp(-Beta_Unadj*(WoSS-Xmin_Unadj));
%%%% Adjusted
WSS = linspace(Xmin_Adj,max(Sample),50);
WCDF = EmpCDF(sum(Sample<Xmin_Adj))*exp(-Beta_Adj*(WSS-Xmin_Adj));
UpperBound = 1-exp(-Beta_Adj.*(max(Sample)-Xmin_Adj));
EmpCDF_Adj = EmpCDF(sum(Sample<Xmin_Adj))*linspace(1,1-UpperBound,sum(Sample>=Xmin_Adj));
EmpCDF_Adj = [EmpCDF(1:sum(Sample<Xmin_Adj)),EmpCDF_Adj];
%%%% Xtick and Ytick
Digi = floor(log(max(Sample))/log(10));
XTick =(floor(min(Sample)/(10^Digi)):ceil(max(Sample)/(10^Digi))).*10^Digi;
YTick = 10.^(floor((log(EmpCDF_Adj(1,end))/log(10))):ceil((log(1)/log(10))));
TextPostionX = linspace(XTick(1,1),XTick(1,end),190);
TextPostionY = 10.^linspace((log(YTick(1,1))/log(10)),(log(YTick(1,end))/log(10)),130);

%% Plotting
figure(1);clf;hold on;
%%%% Unadjusted
subplot(1,2,1);hold on;
plot(Sample,EmpCDF,'.k','markersize',10);
plot(WoSS,WoCDF,'--','color',[1 .0 1],'linewidth',2);
plot(Sample(sum(Sample<Xmin_Unadj)),WoCDF(1,1),'h','color',[1 .0 1],'MarkerFaceColor',[1 .0 1],'markersize',15);
set(gca,'xscale','linear','yscale','log','fontsize',Fontsize);
xlim([XTick(1,1),XTick(1,end)]); ylim([YTick(1,1),YTick(1,end)]);
set(gca,'xtick',XTick,'ytick',YTick);
xlabel(XLabelName,'fontsize',LabelFontsize); ylabel(YLabel,'fontsize',LabelFontsize);
% text(TextPostionX(1,20),TextPostionY(1,25),['P_{KS}: ' sprintf('%0.1f',KS_CSN_Unadj)],'color',[1 .0 1],'fontsize',TextFontsize)
text(TextPostionX(1,20),TextPostionY(1,50),['P_{KS}: ' sprintf('%0.0f',KS_CSN_Unadj) ' / ' sprintf('%0.0f',KS_IF_Unadj)],'color',[1 .0 1],'fontsize',TextFontsize);
text(TextPostionX(1,20),TextPostionY(1,27),['P_{DN}: ' sprintf('%0.0f',DN_CSN_Unadj) ' / ' sprintf('%0.0f',DN_IF_Unadj)],'color',[1 .0 1],'fontsize',TextFontsize);
%%%% Adjusted
subplot(1,2,2);hold on;
plot(Sample,EmpCDF_Adj,'.k','markersize',10);
plot(WoSS,WoCDF,'--','color',[1 .0 1],'linewidth',2);
plot(Sample(sum(Sample<Xmin_Unadj)),WoCDF(1,1),'h','color',[1 .0 1],'MarkerFaceColor',[1 .0 1],'markersize',15);
plot(WSS,WCDF,'--r','linewidth',2);
plot(Sample(sum(Sample<Xmin_Adj)),WCDF(1,1),'h','color','r','MarkerFaceColor','r','markersize',15);
set(gca,'xscale','linear','yscale','log','fontsize',Fontsize);
xlim([XTick(1,1),XTick(1,end)]); ylim([YTick(1,1),YTick(1,end)]);
set(gca,'xtick',XTick,'ytick',YTick);
text(TextPostionX(1,20),TextPostionY(1,50),['P_{KS}: ' sprintf('%0.0f',KS_CSN_Unadj) ' / ' sprintf('%0.0f',KS_IF_Unadj)],'color',[1 .0 1],'fontsize',TextFontsize);
text(TextPostionX(1,20),TextPostionY(1,27),['P_{DN}: ' sprintf('%0.0f',DN_CSN_Unadj) ' / ' sprintf('%0.0f',DN_IF_Unadj)],'color',[1 .0 1],'fontsize',TextFontsize);
text(TextPostionX(1,20),TextPostionY(1,43),['P_{KS}: ' sprintf('%0.0f',KS_CSN_Adj) ' / ' sprintf('%0.0f',KS_IF_Adj)],'color','r','fontsize',TextFontsize);
text(TextPostionX(1,20),TextPostionY(1,20),['P_{DN}: ' sprintf('%0.0f',DN_CSN_Adj) ' / ' sprintf('%0.0f',DN_IF_Adj)],'color','r','fontsize',TextFontsize);
% text(TextPostionX(1,20),TextPostionY(1,18),['P_{KS}: ' sprintf('%0.0f',KS_CSN_Adj) ' / ' sprintf('%0.0f',KS_IF_Adj)],'color','r','fontsize',TextFontsize);
% text(TextPostionX(1,20),TextPostionY(1,10),['P_{DN}: ' sprintf('%0.0f',DN_CSN_Adj) ' / ' sprintf('%0.0f',DN_IF_Adj)],'color','r','fontsize',TextFontsize);
text(TextPostionX(1,end-70),TextPostionY(1,end-10),['\delta: ' sprintf('%0.3f',1-UpperBound)],'color','r','fontsize',TextFontsize);
xlabel(XLabelName,'fontsize',LabelFontsize); ylabel(YLabel,'fontsize',LabelFontsize);


set(gcf,'papersize',[12.5 5])
set(gcf,'paperposition',[.0 .00 12.5 5])
print(gcf,[XLabelName '_PTest'],'-dpdf','-r0')
