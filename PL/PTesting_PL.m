function [Result,PreResult] = PTesting_PL(Sample,SaveDir)
%% This code need to work with QuickKSDNPValueCalculator,PowerLawFittingCompareResult
%%%%%% Sample is empirical data in Nx1 array
Sample = sort(Sample);
Sample = reshape(Sample,[max(size(Sample)),1]);
%% Estimate the parameter X_min and alpah, also determine KS-distance and Distribution noise. 
[WithoutAdjustment,AdjustedForAlpha,AdjustedForAlphaAndKSDist] = PowerLawFittingCompareResult(Sample');
AlphaArraySample = [WithoutAdjustment.PL_Alpha, AdjustedForAlpha.PL_Alpha, AdjustedForAlphaAndKSDist.PL_Alpha];
XminArraySample = [WithoutAdjustment.PL_Xmin, AdjustedForAlpha.PL_Xmin, AdjustedForAlphaAndKSDist.PL_Xmin];
KSDistArraySample = [WithoutAdjustment.PL_KS_Dist, AdjustedForAlpha.PL_KS_Dist, AdjustedForAlphaAndKSDist.PL_KS_Dist];
DistNoiseArraySample = [WithoutAdjustment.PL_DistNoise, AdjustedForAlpha.PL_DistNoise, AdjustedForAlphaAndKSDist.PL_DistNoise];
NumSampleArraySample = [WithoutAdjustment.PL_NumSample, AdjustedForAlpha.PL_NumSample, AdjustedForAlphaAndKSDist.PL_NumSample];

Unadjusted.Parameter.Alpha = AlphaArraySample(1,1);
Unadjusted.Parameter.Xmin = XminArraySample(1,1);
Unadjusted.Parameter.KS = KSDistArraySample(1,1);
Unadjusted.Parameter.DN = DistNoiseArraySample(1,1);
Unadjusted.Parameter.SS = NumSampleArraySample(1,1);

Adjusted.Parameter.Alpha = AlphaArraySample(1,3);
Adjusted.Parameter.Xmin = XminArraySample(1,3);
Adjusted.Parameter.KS = KSDistArraySample(1,3);
Adjusted.Parameter.DN = DistNoiseArraySample(1,3);
Adjusted.Parameter.SS = NumSampleArraySample(1,3);

%% CSN p-testing without adjustment. Resimulate NumSimulation (1000) times of sample using unadjusted parameter. 
NumSimulation = 1000;

XminPWo = XminArraySample(1,1);
AlphaPWo = AlphaArraySample(1,1);
TruncatePWo = Sample(Sample<XminPWo);
SampleSize = size(Sample,1)-size(TruncatePWo,1);

SimulatedSample  = XminPWo*(rand(SampleSize,NumSimulation)).^(1/(1-AlphaPWo));
SimulatedSample = sort([SimulatedSample;repmat(TruncatePWo,[1,NumSimulation])]);

[WithoutAdjustment,AdjustedForAlpha,AdjustedForAlphaAndKSDist] = PowerLawFittingCompareResult(SimulatedSample');
AlphaArrayPWo = [WithoutAdjustment.PL_Alpha, AdjustedForAlpha.PL_Alpha, AdjustedForAlphaAndKSDist.PL_Alpha];
XminArrayPWo = [WithoutAdjustment.PL_Xmin, AdjustedForAlpha.PL_Xmin, AdjustedForAlphaAndKSDist.PL_Xmin];
KSDistArrayPWo = [WithoutAdjustment.PL_KS_Dist, AdjustedForAlpha.PL_KS_Dist, AdjustedForAlphaAndKSDist.PL_KS_Dist];
DistNoiseArrayPWo = [WithoutAdjustment.PL_DistNoise, AdjustedForAlpha.PL_DistNoise, AdjustedForAlphaAndKSDist.PL_DistNoise];
NumSampleArrayPWo = [WithoutAdjustment.PL_NumSample, AdjustedForAlpha.PL_NumSample, AdjustedForAlphaAndKSDist.PL_NumSample];

%% CSN p-testing without adjustment. Resimulate NumSimulation (1000) times of sample using unadjusted parameter. 
XminPW = XminArraySample(1,3);
AlphaPW = AlphaArraySample(1,3);
TruncatePW = Sample(Sample<XminPW);
SampleSize = size(Sample,1)-size(TruncatePW,1);

SimulatedSample  = XminPW*(rand(SampleSize,NumSimulation)).^(1/(1-AlphaPW));
SimulatedSample = sort([SimulatedSample;repmat(TruncatePW,[1,NumSimulation])]);

[WithoutAdjustment,AdjustedForAlpha,AdjustedForAlphaAndKSDist] = PowerLawFittingCompareResult(SimulatedSample');
AlphaArrayPW = [WithoutAdjustment.PL_Alpha, AdjustedForAlpha.PL_Alpha, AdjustedForAlphaAndKSDist.PL_Alpha];
XminArrayPW = [WithoutAdjustment.PL_Xmin, AdjustedForAlpha.PL_Xmin, AdjustedForAlphaAndKSDist.PL_Xmin];
KSDistArrayPW = [WithoutAdjustment.PL_KS_Dist, AdjustedForAlpha.PL_KS_Dist, AdjustedForAlphaAndKSDist.PL_KS_Dist];
DistNoiseArrayPW = [WithoutAdjustment.PL_DistNoise, AdjustedForAlpha.PL_DistNoise, AdjustedForAlphaAndKSDist.PL_DistNoise];
NumSampleArrayPW = [WithoutAdjustment.PL_NumSample, AdjustedForAlpha.PL_NumSample, AdjustedForAlphaAndKSDist.PL_NumSample];

%% Calculate p-value
KS_SSDependencePower = 0.4924;
DN_SSDependencePower = 0.4946;
%%% CSN unadjusted
KS_Sample_Unadj = KSDistArraySample(1,1);
DN_Sample_Unadj = DistNoiseArraySample(1,1)-1/sqrt(2);
SS_Sample_Unadj = NumSampleArraySample(1,1);
KS_CSN_Unadj = KSDistArrayPWo(:,1);
DN_CSN_Unadj = DistNoiseArrayPWo(:,1)-1/sqrt(2);
SS_CSN_Unadj = NumSampleArrayPWo(:,1);

KSPValue_CSN_Unadj = 100-100*sum(KS_CSN_Unadj.*(SS_CSN_Unadj.^KS_SSDependencePower)<KS_Sample_Unadj.*(SS_Sample_Unadj.^KS_SSDependencePower))/NumSimulation;
DNPValue_CSN_Unadj = 100-100*sum(DN_CSN_Unadj.*(SS_CSN_Unadj.^DN_SSDependencePower)<DN_Sample_Unadj.*(SS_Sample_Unadj.^DN_SSDependencePower))/NumSimulation;

Unadjusted.PValue.CSN_KS = KSPValue_CSN_Unadj;
Unadjusted.PValue.CSN_DN = DNPValue_CSN_Unadj;
%%% Inversion formulae unadjusted
[PKS,PDN] = QuickKSDNPValueCalculator(KS_Sample_Unadj,DN_Sample_Unadj+1/sqrt(2),SS_Sample_Unadj);
Unadjusted.PValue.IF_KS = 100-PKS;
Unadjusted.PValue.IF_DN = 100-PDN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CSN adjusted
KS_Sample_adj = KSDistArraySample(1,3);
DN_Sample_adj = DistNoiseArraySample(1,3)-1/sqrt(2);
SS_Sample_adj = NumSampleArraySample(1,3);
KS_CSN_adj = KSDistArrayPW(:,3);
DN_CSN_adj = DistNoiseArrayPW(:,3)-1/sqrt(2);
SS_CSN_adj = NumSampleArrayPW(:,3);

KSPValue_CSN_adj = 100-100*sum(KS_CSN_adj.*(SS_CSN_adj.^KS_SSDependencePower)<KS_Sample_adj.*(SS_Sample_adj.^KS_SSDependencePower))/NumSimulation;
DNPValue_CSN_adj = 100-100*sum(DN_CSN_adj.*(SS_CSN_adj.^DN_SSDependencePower)<DN_Sample_adj.*(SS_Sample_adj.^DN_SSDependencePower))/NumSimulation;

Adjusted.PValue.CSN_KS = KSPValue_CSN_adj;
Adjusted.PValue.CSN_DN = DNPValue_CSN_adj;
%%% Inversion formulae adjusted
[PKS,PDN] = QuickKSDNPValueCalculator(KS_Sample_adj,DN_Sample_adj+1/sqrt(2),SS_Sample_adj);
Adjusted.PValue.IF_KS = 100-PKS;
Adjusted.PValue.IF_DN = 100-PDN;

Result.Sample = Sample;
Result.Adjusted = Adjusted;
Result.Unadjusted = Unadjusted;
%% Organized Data with estimated parameters details
PreResult.Sample = Sample;
PreResult.AlphaArraySample = AlphaArraySample;
PreResult.XminArraySample = XminArraySample;
PreResult.KSDistArraySample = KSDistArraySample;
PreResult.DistNoiseArraySample = DistNoiseArraySample;
PreResult.NumSampleArraySample = NumSampleArraySample;

PreResult.AlphaArrayPWo = AlphaArrayPWo;
PreResult.XminArrayPWo = XminArrayPWo;
PreResult.KSDistArrayPWo = KSDistArrayPWo;
PreResult.DistNoiseArrayPWo = DistNoiseArrayPWo;
PreResult.NumSampleArrayPWo = NumSampleArrayPWo;

PreResult.AlphaArrayPW = AlphaArrayPW;
PreResult.XminArrayPW = XminArrayPW;
PreResult.KSDistArrayPW = KSDistArrayPW;
PreResult.DistNoiseArrayPW = DistNoiseArrayPW;
PreResult.NumSampleArrayPW = NumSampleArrayPW;

save([SaveDir '/Result_CSNTest.mat'],'Result','PreResult');