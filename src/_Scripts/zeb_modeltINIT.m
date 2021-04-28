%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% METADATA %%%%
%
% zeb_modeltINIT.m
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
% [1] ZebraGEM/Data/ZebraGEM2.1.xlsx
%     Generalised GEM.
%
% [2] zeb/Results/DESeq2/all/
%     DESeq2 results.
%
% [3] zeb/Results/ConsensusExpression/
%     Mean expression values for clusters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%
%
% [1] zeb/Results/GEM/ZebraGEM2.1-brain-XXX.sbml.xml
% [2] zeb/Results/GEM/ZebraGEM2.1-brain-XXX.xlsx
%     Context-specific GEMs.
%
% [3] zeb/Results/GEM/fluxes_WWW_XXX_vs_YYY.txt
%     Flux comparison between genotypes X and Y in tissue W.
%
% [4] zeb/Results/GEM/repMets_WWW_XXX_vs_YYY.tsv
%     Reporter metabolite analysis results between genotypes X and Y in tissue W.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%
%% Initialise

% Import reference model
zebmodel = importExcelModel('ZebraGEM/Data/ZebraGEM2.1.xlsx', false); % corrected model
tasksLC = parseTaskList('GEM/Data/common_tasks_growth_RPMI1640.xlsx');     % Load in the task list

% Set objective function to ATP synthesis
zebmodel = setParam(zebmodel, 'obj', 'BIO_L', 0);
zebmodel = setParam(zebmodel, 'obj', 'ATPS4m', 1);

params.TimeLimit = 30000;                                                                       % Gurobi: time limit in sec
params.MIPGap = 0.02;                                                                           % Gurobi: relative MIP optimality gap

% I/O
nClust = 15;                                                                                       % Number of models (clsuters) in the set
mod(1).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_full_dko.tsv';
mod(2).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_full_het.tsv';
mod(3).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_full_wt.tsv';
mod(4).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_brain_dko.tsv';
mod(5).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_brain_het.tsv';
mod(6).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_brain_wt.tsv';
mod(7).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_muscle_dko.tsv';
mod(8).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_muscle_het.tsv';
mod(9).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_muscle_wt.tsv';
mod(10).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_skin_dko.tsv';
mod(11).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_skin_het.tsv';
mod(12).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_skin_wt.tsv';
mod(13).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_liver_dko.tsv';
mod(14).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_liver_het.tsv';
mod(15).infile = '/Users/aceid/PhD/2Year2/2ZebTel/MeanCts/Results/meanCts_liver_wt.tsv';
mod(1).name = 'Zebrafish-Full-dko';
mod(2).name = 'Zebrafish-Full-het';
mod(3).name = 'Zebrafish-Full-wt';
mod(4).name = 'Zebrafish-Brain-dko';
mod(5).name = 'Zebrafish-Brain-het';
mod(6).name = 'Zebrafish-Brain-wt';
mod(7).name = 'Zebrafish-Muscle-dko';
mod(8).name = 'Zebrafish-Muscle-het';
mod(9).name = 'Zebrafish-Muscle-wt';
mod(10).name = 'Zebrafish-Skin-dko';
mod(11).name = 'Zebrafish-Skin-het';
mod(12).name = 'Zebrafish-Skin-wt';
mod(13).name = 'Zebrafish-Liver-dko';
mod(14).name = 'Zebrafish-Liver-het';
mod(15).name = 'Zebrafish-Liver-wt';
outfile_dir = 'zeb/Results/GEM/';
infile_DESeq2 = 'zeb/Results/DESeq2/all/sep/entrez/';

%% tINIT by cluster

for i = 1:nClust

    % Import cluster-wise expression data
    mod(i).consExpr = readTXT(mod(i).infile);
    mod(i).genes = mod(i).consExpr(2:end,1);
    mod(i).values = str2double(mod(i).consExpr(2:end,2));

    % Prepare tINIT input
    mod(i).tINITinput = exp2tINITinput(mod(i).genes,mod(i).values,mod(i).name);

    % Import required task list

    try
        % Do tINIT
        mod(i).model = getINITModel(zebmodel, mod(i).name, mod(i).name, mod(i).tINITinput, [], [], [], [], true, tasksLC, params,[]);

        % Export the model
        exportModel(mod(i).model, strcat(outfile_dir, mod(i).name, '.sbml.xml'));
        exportToExcelFormat(mod(i).model, strcat(outfile_dir, mod(i).name, '.xlsx'));
    catch
        warning(['An error occurred. Infeasible model?']);               %  if an error is encountered, don't stop execution
    end

end

%% Simulate the models

for i = 1:nClust

    % Constraints
    mod(i).model = mod(i).modeltemp;
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_gln_L_e',0.013,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_gln_L_e',-0.025,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_acac_e',0.012,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_acac_e',0.0015,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'r1391',100,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'r1391',0,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'G3PDm',100,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'G3PDm',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_h2o2_e',100,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_h2o2_e',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_hco3_e',100,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_hco3_e',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_lac_L_e',0.0058,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_lac_L_e',-0.079,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_adrnl_e',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_adrnl_e',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_dopa_e',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_dopa_e',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_nrpphr_e',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_nrpphr_e',0,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_ile_e',0.0041,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_ile_e',-0.0004,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_hxan_e',0.00045,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_hxan_e',-0.00045,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ins_e',0.00045,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ins_e',-0.00045,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ser_L_e',0.011,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ser_L_e',-0.0016,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_val_e',0.011,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_val_e',-0.005,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_asn_L_e',0.0009,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_asn_L_e',-0.0037,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_pro_L_e',0.0079,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_pro_L_e',-0.0066,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_pyr_e',0.0058,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_pyr_e',-0.007,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ala_L_e',100,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ala_L_e',-0.0079,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_gly_e',0.0053,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_gly_e',-0.0086,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_lys_e',0.0005,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_lys_e',-0.011,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_glu_e',0.0044,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_glu_e',-0.0047,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_glc_bD_e',0.29,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_glc_bD_e',0.196,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_phe_e',0,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_phe_e',-100,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_trp_e',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_trp_e',-100,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'GLCt1r',0.19,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'GLCt1r',0.16,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_thr_e',0,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_thr_e',-0.0008,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_met__L_e',0,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_met__L_e',-0.0017,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_tyr_L_e',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_tyr_L_e',-0.0017,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_arg_e',0.004,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_arg_e',0,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_his__L_e',0,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_his__L_e',-0.0025,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_orn_e',0.0048,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_orn_e',-0.0041,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_leu_e',0.0062,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_leu_e',-0.0011,'u');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_o2_e',2.256,'l');
%     mod(i).model = changeRxnBounds(mod(i).model, 'EX_o2_e',1.351,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_co2_e',-0.515,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_co2_e',-0.53,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_3aib_e',0.00023,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_4abut_e',0.0018,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_4abut_e',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ac_e',0.0013,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ac_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_acald_e',0.0014,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_acald_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ala_B_e',1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_ala_B_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_bhb_e',0.016,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_bhb_e',0.001,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_btn_e',1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_btn_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_crn_e',1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_crn_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_cys_L_e',0.0086,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_cys_L_e',-0.0033,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_duri_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_duri_e',1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_etoh_e',0,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_etoh_e',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_fe2_e',1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_fe2_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_nh3_e',1000,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_nh3_e',0,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_pe_hs_e',0.00005,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'EX_pe_hs_e',-1000,'u');
    mod(i).model = changeRxnBounds(mod(i).model, 'ps_hs_mem',0.00005,'l');
    mod(i).model = changeRxnBounds(mod(i).model, 'ps_hs_mem',-1000,'u');

    % pFBA
    %mod(i).model = simplifyModel(mod(i).model, true, false, true, true);                        % Simplify the model
    mod(i).model = setParam(mod(i).model, 'obj', 'ATPS4m', 1);                                % Set objective function to be ATP synthesis
    mod(i).sol = solveLP(mod(i).model);                                                   % pFBA

end

%% Compare results between clusters

% pFBA
diary 'zeb/Results/GEM/fluxes_full_dko_vs_wt.txt'
followChanged(zebmodel, mod(1).sol.x, mod(3).sol.x);
diary 'zeb/Results/GEM/fluxes_full_het_vs_wt.txt'
followChanged(zebmodel, mod(2).sol.x, mod(3).sol.x);
diary 'zeb/Results/GEM/fluxes_brain_dko_vs_wt.txt'
followChanged(zebmodel, mod(4).sol.x, mod(6).sol.x);
diary 'zeb/Results/GEM/fluxes_brain_het_vs_wt.txt'
followChanged(zebmodel, mod(5).sol.x, mod(6).sol.x);
diary 'zeb/Results/GEM/fluxes_muscle_wt_vs_dko.txt'
followChanged(zebmodel, mod(9).sol.x, mod(7).sol.x);
diary 'zeb/Results/GEM/fluxes_muscle_wt_vs_het.txt'
followChanged(zebmodel, mod(9).sol.x, mod(8).sol.x);
diary 'zeb/Results/GEM/fluxes_skin_dko_vs_wt.txt'
followChanged(zebmodel, mod(10).sol.x, mod(12).sol.x);
diary 'zeb/Results/GEM/fluxes_skin_het_vs_wt.txt'
followChanged(zebmodel, mod(11).sol.x, mod(12).sol.x);
diary 'zeb/Results/GEM/fluxes_liver_dko_vs_wt.txt'
followChanged(zebmodel, mod(13).sol.x, mod(15).sol.x);
diary 'zeb/Results/GEM/fluxes_liver_wt_vs_het.txt'
followChanged(zebmodel, mod(15).sol.x, mod(14).sol.x);
diary off

% Reporter metabolite analysis

%
repMet(1).comparison = 'Full_dko_vs_wt';
repMet(1).infile = strcat(infile_DESeq2, 'Full_dko_vs_wt_entrez.tsv');
repMet(1).outfile = strcat(outfile_dir, 'repMets_Full_dko_vs_wt.txt');
repMet(2).comparison = 'Full_het_vs_wt';
repMet(2).infile = strcat(infile_DESeq2, 'Full_het_vs_wt_entrez.tsv');
repMet(2).outfile = strcat(outfile_dir, 'repMets_Full_het_vs_wt.txt');
repMet(3).comparison = 'Full_dko_vs_het';
repMet(3).infile = strcat(infile_DESeq2, 'Full_dko_vs_het_entrez.tsv');
repMet(3).outfile = strcat(outfile_dir, 'repMets_Full_dko_vs_het.txt');
repMet(4).comparison = 'Brain_dko_vs_wt';
repMet(4).infile = strcat(infile_DESeq2, 'Brain_dko_vs_wt_entrez.tsv');
repMet(4).outfile = strcat(outfile_dir, 'repMets_Brain_dko_vs_wt.txt');
repMet(5).comparison = 'Brain_het_vs_wt';
repMet(5).infile = strcat(infile_DESeq2, 'Brain_het_vs_wt_entrez.tsv');
repMet(5).outfile = strcat(outfile_dir, 'repMets_Brain_het_vs_wt.txt');
repMet(6).comparison = 'Brain_dko_vs_het';
repMet(6).infile = strcat(infile_DESeq2, 'Brain_dko_vs_het_entrez.tsv');
repMet(6).outfile = strcat(outfile_dir, 'repMets_Brain_dko_vs_het.txt');
repMet(7).comparison = 'Muscle_dko_vs_wt';
repMet(7).infile = strcat(infile_DESeq2, 'Muscle_dko_vs_wt_entrez.tsv');
repMet(7).outfile = strcat(outfile_dir, 'repMets_Muscle_dko_vs_wt.txt');
repMet(8).comparison = 'Muscle_het_vs_wt';
repMet(8).infile = strcat(infile_DESeq2, 'Muscle_het_vs_wt_entrez.tsv');
repMet(8).outfile = strcat(outfile_dir, 'repMets_Muscle_het_vs_wt.txt');
repMet(9).comparison = 'Muscle_dko_vs_het';
repMet(9).infile = strcat(infile_DESeq2, 'Muscle_dko_vs_het_entrez.tsv');
repMet(9).outfile = strcat(outfile_dir, 'repMets_Muscle_dko_vs_het.txt');
repMet(10).comparison = 'Skin_dko_vs_wt';
repMet(10).infile = strcat(infile_DESeq2, 'Skin_dko_vs_wt_entrez.tsv');
repMet(10).outfile = strcat(outfile_dir, 'repMets_Skin_dko_vs_wt.txt');
repMet(11).comparison = 'Skin_het_vs_wt';
repMet(11).infile = strcat(infile_DESeq2, 'Skin_het_vs_wt_entrez.tsv');
repMet(11).outfile = strcat(outfile_dir, 'repMets_Skin_het_vs_wt.txt');
repMet(12).comparison = 'Skin_dko_vs_het';
repMet(12).infile = strcat(infile_DESeq2, 'Skin_dko_vs_het_entrez.tsv');
repMet(12).outfile = strcat(outfile_dir, 'repMets_Skin_dko_vs_het.txt');
repMet(13).comparison = 'Liver_dko_vs_wt';
repMet(13).infile = strcat(infile_DESeq2, 'Liver_dko_vs_wt_entrez.tsv');
repMet(13).outfile = strcat(outfile_dir, 'repMets_Liver_dko_vs_wt.txt');
repMet(14).comparison = 'Liver_het_vs_wt';
repMet(14).infile = strcat(infile_DESeq2, 'Liver_het_vs_wt_entrez.tsv');
repMet(14).outfile = strcat(outfile_dir, 'repMets_Liver_het_vs_wt.txt');
repMet(15).comparison = 'Liver_dko_vs_het';
repMet(15).infile = strcat(infile_DESeq2, 'Liver_dko_vs_het_entrez.tsv');
repMet(15).outfile = strcat(outfile_dir, 'repMets_Liver_dko_vs_het.txt');

for i = 1:nClust
    repMet(i).DESeq2 = readTXT(repMet(i).infile);
    repMet(i).genes = repMet(i).DESeq2(2:end,1);
    repMet(i).l2FC = str2double(repMet(i).DESeq2(2:end,3));
    repMet(i).pvalues = str2double(repMet(i).DESeq2(2:end,6));
    repMet(i).repMets = reporterMetabolites(zebmodel, repMet(i).genes, repMet(i).pvalues, true, repMet(i).outfile, repMet(i).l2FC);
end
