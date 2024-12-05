clc, clearvars;

% Adding nanoCOBRA to the path

addpath(genpath('nanoCobratoolbox'))

% Reading the Oryza sativa model

model = readCbModel('o_sativa_core.mat');

% Simulating seed cells in aerobic and anaerobic conditions

model_copy = changeObjective(model,'Coleoptile_Biomass');
model_copy = changeRxnBounds(model_copy,{'Ex_Sucrose'},-1,'b');

o2_rates = {-1, -2, -3, -3.35, -4};

aerobic_fluxes = table(model_copy.rxns,'VariableNames', {'Rxn'});

for rate = 1:length(o2_rates)

    temp_model = changeRxnBounds(model_copy,{'Ex_O2'},o2_rates{rate},'l');
    temp_sol = optimizeCbModel(temp_model);
    aerobic_fluxes.(num2str(-o2_rates{rate})) = temp_sol.v;

end

% Simulating seed cells' sucrose batch cultures in aerobic and anaerobic conditions

model_copy = changeObjective(model,'Coleoptile_Biomass');
o2_rates = {-3.35, 0};
suc_rates = {-8.224, -9.5};
glc_rates = {3.472, 0};
fru_rates = {4.167, 5.556};
experiments = {'Sucrose Aerobic', 'Sucrose Anaerobic'};

sugar_fluxes = table(model_copy.rxns,'VariableNames', {'Rxn'});

for rate = 1:length(o2_rates)

    temp_model = changeRxnBounds(model_copy,{'Ex_O2'},o2_rates{rate},'l');
    temp_model = changeRxnBounds(temp_model,{'Ex_Sucrose'},suc_rates{rate},'b');
    temp_model = changeRxnBounds(temp_model,{'Ex_Glucose'},glc_rates{rate},'b');
    temp_model = changeRxnBounds(temp_model,{'Ex_Fructose'},fru_rates{rate},'b');
    temp_sol = optimizeCbModel(temp_model);
    sugar_fluxes.(experiments{rate}) = temp_sol.v;

end

% Simulating seed cells' glucose batch cultures in aerobic and anaerobic conditions

model_copy = changeObjective(model,'Coleoptile_Biomass');
o2_rates = {-3.35, 0};
glc_rates = {-2.778, -19.444};
experiments = {'Glucose Aerobic', 'Glucose Anaerobic'};

for rate = 1:length(o2_rates)

    temp_model = changeRxnBounds(model_copy,{'Ex_O2'},o2_rates{rate},'l');
    temp_model = changeRxnBounds(temp_model,{'Ex_Glucose'},glc_rates{rate},'b');
    temp_sol = optimizeCbModel(temp_model);
    sugar_fluxes.(experiments{rate}) = temp_sol.v;

end

% Simulating leaf cells in varying O2 and CO2 conditions

model_copy = changeObjective(model,'Straw_Biomass');
model_copy = changeRxnBounds(model_copy,{'Ex_photon'},-100,'b');
o2_rates = {0,-1,-1,-1,-1,-1,-1};
co2_rates = {-1,-1,-2,-4,-6,-8,-10};
experiments = {'1:0','1:1','2:1','4:1','6:1','8:1','10:1'};

photorespire_fluxes = table(model_copy.rxns,'VariableNames', {'Rxn'});

for rate = 1:length(o2_rates)
    temp_model = changeRxnBounds(model_copy,{'Ex_O2'},-o2_rates{rate},'u');
    temp_model = changeRxnBounds(temp_model,{'Ex_CO2'},-co2_rates{rate},'u');
    temp_sol = optimizeCbModel(temp_model);
    photorespire_fluxes.(experiments{rate}) = temp_sol.v;

end

% Flux variability analysis

o2_rates = {-3.35,0,-3.35};
growth_rates = {0.3750,0.1133,0.174};
objective = {'Coleoptile_Biomass','Coleoptile_Biomass','Straw_Biomass'};
experiments = {'Aerobic Seed','Anaerobic Seed','Aerobic Leaf'};

flux_variability = table(model_copy.rxns,'VariableNames', {'Rxn'});

for rate = 1:length(o2_rates)

    model_copy = changeObjective(model,objective{rate});
    model_copy = changeRxnBounds(model_copy,{'Ex_O2'},o2_rates{rate},'b');
    model_copy = changeRxnBounds(model_copy,{objective{rate}},growth_rates{rate},'b');
    [minFlux, maxFlux] = fluxVariability(model_copy,90);
    flux_variability.(strcat(experiments{rate}," min")) = minFlux; 
    flux_variability.(strcat(experiments{rate}," max")) = maxFlux; 

end

filename = 'be21b004_Results.xlsx';
writetable(aerobic_fluxes,filename,'Sheet',"seed-aerobic-anaerobic",'WriteMode','overwritesheet')
writetable(sugar_fluxes,filename,'Sheet',"seed-sucrose-glucose",'WriteMode','overwritesheet')
writetable(photorespire_fluxes,filename,'Sheet',"leaf-photorespiration",'WriteMode','overwritesheet')
writetable(flux_variability,filename,'Sheet',"fva",'WriteMode','overwritesheet')
