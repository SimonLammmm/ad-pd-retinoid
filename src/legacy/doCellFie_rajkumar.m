%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% METADATA %%%%
%
% doCellFie_rajkumar.m
% Version: 1.0.3 (2020-10-20)
%
% Author: Simon Lam
% Institution: King's College London
% Contact: simon.1.lam@kcl.ac.uk
%
% Description:
% Conduct cellular functions inference.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEGAL %%%%
%
% Copyright © 2020 King's College London
%
% This work is licensed under the Creative Commons Attribution 4.0 International Licence. To view
% a copy of this license,  visit http://creativecommons.org/licences/by/4.0/  or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
%
% Permission is hereby granted,  free of charge, to any person  obtaining a copy of this software
% and  associated  documentation  files  (the  "Software"),  to  deal  in  the  Software  without
% restriction,  including without  limitation the  rights to use,  copy, modify,  merge, publish,
% distribute, and/or sell copies  of the Software, and to permit persons  to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above  copyright notice  and this  permission notice  shall be included  with all copies or
% substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED  "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS OR IMPLIED, INCLUDING
% BUT NOT  LIMITED TO THE WARRANTIES  OF MERCHANTABILITY,  FITNESS FOR A  PARTICULAR PURPOSE, AND
% NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS  OR COPYRIGHT HOLDERS BE  LIABLE FOR ANY CLAIM,
% DAMAGES, OR  OTHER LIABILITY,  WHETHER IN AN  ACTION OF CONTRACT,  TORT, OR  OTHERWISE, ARISING
% FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF OR DEALING IN THE SOFTWARE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%
%
% [1] supplementary_rajkumar/Results/ConsensusExpression/countsConsensusEntrez.tsv
%     Mean expression values per cluster (Entrez format)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%
%
% [1] supplementary_rajkumar/Results/CellFie/CellFie.csv
%     CellFie results.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%
%%
THIS_DIR = pwd;
initCellFie;
cd(THIS_DIR);

%% Load dataset
inpath = 'supplementary_rajkumar/Results/ConsensusExpression/countsConsensusEntrez.tsv';
raw = readTXT(inpath);
data.Tissue = raw(1,2:end);
data.gene = raw(2:end,1);
data.value = str2double(raw(2:end,2:end));
outpath = 'supplementary_rajkumar/Results/CellFie/';

mkdir(outpath);

%% Define the number of samples (equal to the column number of the expression matrix)
SampleNumber = length(data.Tissue);

%% Define the reference genome-scale models you want to use (all listed in the test/suite)
ref='MT_recon_2_2_entrez.mat';

%% Define the type parameters of the method
param.ThreshType='local';
param.LocalThresholdType='minmaxmean';
param.percentile_or_value='value';
param.value_low=25;
param.value_high=75;

%% Run CellFie
[out.score, out.score_binary, out.taskInfos, out.detailScoring] = CellFie(data,SampleNumber,ref,param);

%% Parse output
out2.essentialReaction = out.detailScoring(:,5);
out2.expressionScore = out.detailScoring(:,3:8:SampleNumber * 8);
out2.expressionScore = num2cell(cellfun(@abs, out2.expressionScore));
out2.Tissue = data.Tissue;

%% Save output
out3 = cell2table([out2.essentialReaction out2.expressionScore], 'VariableNames', raw(1,1:end));
writetable(out3, strcat(outpath,'CellFie.csv'));
