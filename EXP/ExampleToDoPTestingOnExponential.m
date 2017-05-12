%%%%% This code need to work with QuickKSDNPValueCalculator,PowerLawFittingCompareResult,Plot_PL
%%%%%  The code will take hours or days to finish
%% Perform P-testing
load('EXP_Sample.mat');
Sample = EXP_Sample;
[Result,PreResult] = PTesting_EXP(Sample,cd);
%% Plot Result
load('Result_CSNTest.mat');
Plot_EXP(Result,'EXP Sample');
