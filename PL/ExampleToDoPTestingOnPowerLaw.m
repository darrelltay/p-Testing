%%%%% This code need to work with QuickKSDNPValueCalculator,PowerLawFittingCompareResult,Plot_PL
%%%%%  The code will take hours or days to finish
%% Perform P-testing
% load('PL_Sample.mat');
% Sample = PL_Sample;
% [Result,PreResult] = PTesting_PL(Sample,cd);
%% Plot Result
load('Result_CSNTest.mat');
Plot_PL(Result,'PL Sample');
