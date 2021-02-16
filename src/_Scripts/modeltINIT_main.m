%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% METADATA %%%%
%
% modeltINIT_main.m
% Version: 1.1.0 (2020-11-02)
%
% Author: Simon Lam
% Institution: King's College London
% Contact: simon.1.lam@kcl.ac.uk
%
% Description:
% Generate context-specific GEMs and analyse them.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE REQUIREMENTS %%%%
%
% [1] GEM/Data/common_tasks_growth_RPMI1640.xlsx
%     Common growth task list.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%
%
% [1] iBrain2845/Results/iBrain2845.sbml.xlsx
%     Generalised GEM.
%
% [2] main/Results/DESeq2/
%     DESeq2 results.
%
% [3] main/Results/ConsensusExpression/
%     Mean expression values for clusters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%
%
% [1] main/Results/GEM/iADPDX.sbml.xml
% [2] main/Results/GEM/iADPDX.xlsx
%     Context-specific GEMs.
%
% [3] main/Results/GEM/fluxes_cX_vs_cY.txt
%     Flux comparison between clusters X and Y.
%
% [4] main/Results/GEM/repMets_cX_vs_cY.tsv
%     Reporter metabolite analysis results between clusters X and Y.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%

% Import
refmodel = importExcelModel('iBrain2845/Results/iBrain2845.sbml.xlsx', false);
refmodel.rxnNames = refmodel.rxns;
tasksLC = parseTaskList('GEM/Data/common_tasks_growth_RPMI1640.xlsx');         % Load in the task list
infile_DESeq2 = 'main/Results/DESeq2/';

setRavenSolver('gurobi')
params.TimeLimit = 30000;                                                                       % Gurobi: time limit in sec
params.MIPGap = 0.02;                                                                           % Gurobi: relative MIP optimality gap

% I/O
nClust = 4;                                                                                       % Number of models (clusters) in the set
infile_dir = 'main/Results/ConsensusExpression/';
outfile_dir = 'main/Results/GEM/';

%% tINIT by cluster

for i = 1:nClust

    % Import cluster-wise expression data
    mod(i).consExpr = readTXT(strcat(infile_dir, 'cluster_', num2str(i), '_consensus_expression.tsv'));
    mod(i).genes = mod(i).consExpr(2:end,1);
    mod(i).values = str2double(mod(i).consExpr(2:end,2));

    % Prepare tINIT input
    mod(i).modelName = strcat('Cluster', num2str(i));
    mod(i).tINITinput = exp2tINITinput(mod(i).genes,mod(i).values,mod(i).modelName);

    %try
        % Do tINIT
        mod(i).model = getINITModel(refmodel, mod(i).modelName, mod(i).modelName, mod(i).tINITinput, [], [], [], [], true, tasksLC, params,[]);

        % Export the model
        exportModel(mod(i).model, strcat(outfile_dir, 'iADPD', num2str(i), '.sbml.xml'));
        exportToExcelFormat(mod(i).model, strcat(outfile_dir, 'iADPD', num2str(i), '.xlsx'));
    %catch
    %    warning(['An error occurred. Infeasible model?']);               %  if an error is encountered, don't stop execution
    %end

end

%% Simulate the models

for i = 1:nClust

    % Constraints
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9063',0.013,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9063',-0.025,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9132',0.012,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9132',0.0015,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_5399',100,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_5399',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_5400',100,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_5400',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_0482',100,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_0482',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_h2o2_s',100,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_h2o2_s',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9078',100,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9078',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9135',0.0058,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9135',-0.079,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9095',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9095',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9092',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9092',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9093',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9093',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9039',0.0041,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9039',-0.0004,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_hxan_s',0.00045,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_hxan_s',-0.00045,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9361',0.00045,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9361',-0.00045,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9069',0.011,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9069',-0.0016,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9046',0.011,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9046',-0.005,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9062',0.0009,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9062',-0.0037,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9068',0.0079,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9068',-0.0066,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9133',0.0058,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9133',-0.007,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9061',100,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9061',-0.0079,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9067',0.0053,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9067',-0.0086,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9041',0.0005,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9041',-0.011,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9071',0.0044,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9071',-0.0047,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9034',0.29,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9034',0.196,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9043',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9043',-100,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9045',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9045',-100,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_5029',0.19,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_5029',0.16,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9044',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9044',-0.0008,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9042',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9042',-0.0017,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9064',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9064',-0.0017,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9066',0.004,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9066',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9038',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9038',-0.0025,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9087',0.0048,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9087',-0.0041,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9040',0.0062,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9040',-0.0011,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9048',2.256,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9048',1.351,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9058',-0.515,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9058',-0.53,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_3aib_s',0.00023,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9091',0.0018,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9091',0,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9086',0.0013,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9086',-1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_acald_s',0.0014,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_acald_s',-1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9260',1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9260',-1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9134',0.016,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9134',0.001,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_crn_s',1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_crn_s',-1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9065',0.0086,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9065',-0.0033,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_duri_s',-1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_duri_s',1000,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9099',0,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9099',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9076',1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9076',-1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9073',1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'HMR_9073',0,'l');

    % pFBA
    mod(i).model = simplifyModel(mod(i).model, true, false, true, true);                        % Simplify the model
    mod(i).model = setParam(mod(i).model, 'obj', 'HMR_6916', 1);                                % Set objective function to be ATP synthesis
    mod(i).sol = solveLP(mod(i).model, true);                                                   % pFBA

end

%% Compare results between clusters

% pFBA
diary 'main/Results/GEM/fluxes_c4_vs_c1.txt'
followChanged(refmodel, mod(4).sol.x, mod(1).sol.x)
diary 'main/Results/GEM/fluxes_c4_vs_c2.txt'
followChanged(refmodel, mod(4).sol.x, mod(2).sol.x)
diary 'main/Results/GEM/fluxes_c3_vs_c4.txt'
followChanged(refmodel, mod(3).sol.x, mod(4).sol.x)
diary off

% Reporter metabolite analysis
%
for i = 1:nClust-1
    repMet.infile = strcat(infile_DESeq2, 'allDEG_c', num2str(i), '_vs_c', num2str(nClust), '.tsv');
    repMet.DESeq2 = readTXT(repMet.infile);
    repMet.genes = repMet.DESeq2(2:end,1);
    repMet.l2FC = str2double(repMet.DESeq2(2:end,3));
    repMet.pvalues = str2double(repMet.DESeq2(2:end,6));
    repMet.outfile = strcat(outfile_dir, 'repMets_c', num2str(i), '_vs_c', num2str(nClust), '.txt');
    repMet.repMets = reporterMetabolites(refmodel, repMet.genes, repMet.pvalues, true, repMet.outfile, repMet.l2FC);
end
