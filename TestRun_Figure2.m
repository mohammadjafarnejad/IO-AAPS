%% Produce Figure 2 using the provided model as a sample case 
%  This is a selected partial response case from the parameter sensitivity 
%  analysis cases. This case is different from the baseline case that 
%  contains the baseline parameters provided in the supplements in SBML
%  format and listed in an excel file.
%
%
% Created: April 02, 2019 (Mohammad Jafarnejad)
% Last Modified: April 02, 2019 (MJ)

%% Clear the workspace and close all the figures
clear all
close all
clc

%% Load the model
load('Jafarnejad et al AAPS 2019.mat')

%% Add cell unit used by the model
unitObj = sbiounit('cell', 'molecule', 1);
sbioaddtolibrary(unitObj);

%% Get the list of species, parameters and reactions
for i = 1:size(model.Parameters,1)
    Parameters(i,1) = {model.Parameters(i).Name};
    Parameters(i,2) = {model.Parameters(i).Value};
    Parameters(i,3) = {model.Parameters(i).Tag};
end
% ------------------------------------------------------------------- 
for i = 1:size(model.Species,1)
    Species(i,1) = {model.Species(i).Name};
    Species(i,2) = {model.Species(i).InitialAmount};
    Species(i,3) = {model.Species(i).Parent.Name};
    Species(i,4) = {model.Species(i).Tag};
end
% ------------------------------------------------------------------- 
for i = 1:size(model.reaction,1)
    ReactionRates(i,1) = {model.reaction(i).Reaction};
    ReactionRates(i,2) = {model.reaction(i).ReactionRate};
    ReactionRates(i,3) = {model.reaction(i).Tag};
end

%% Run to reach steady state
model.Reactions(find(strcmp(ReactionRates(:,3),'SteadyState'))).Active = 1;
model.Parameters(find(strcmp(Parameters(:,1),'dose_mg_Nivo'  ))).Value = 0;
simDataSS  = sbiosimulate(model);

%% Replace the species with steady state values
    speciesSS = simDataSS.data(end,:)'; 
    for k = 1:size(model.Species,1)
        model.Species(k).InitialAmount = speciesSS(k);
    end
    
%% Revert the model that contains initial condition back to dosing mode and turn on dosing at 30 days
    model.Reactions(find(strcmp(ReactionRates(:,3),'SteadyState'))).Active = 0;
    model.Parameters(find(strcmp(Parameters(:,1),'dose_mg_Nivo'  ))).Value = 3;
    model.Parameters(find(strcmp(Parameters(:,1),'t_init_Nivo'   ))).Value = 30;
    simData  = sbiosimulate(model);
    
%% Plot the subplots in the figure
f = figure; 
set(f,'Position', [600 50 600 900]); 
linesize = 2; fontsize = 12;

subplot(3,2,1); hold on; box on; clear i; 
i =[find( strcmp(simData.DataNames(:,1),'D_tum_app') )]; 
plot(simData.time/30, simData.data(:,i), 'LineWidth', linesize)
ylabel('Tumor Diameter (mm)','Fontsize',fontsize); 
xlabel('Time (month)','Fontsize',fontsize);
xlim([min(simData.time/30) max(simData.time/30)]);
ylim([0 100])
set(gca,'TickLength',[0.03, 0.015])
set(gca,'XTick',[0, 6, 12])

subplot(3,2,2); hold on; box on; clear i; 
i =[find(strcmp(Species(:,1),'APC') & strcmp(Species(:,3),'Tum'));...
    find(strcmp(Species(:,1),'mAPC')& strcmp(Species(:,3),'Tum'));...
    find(strcmp(Species(:,1),'APC') & strcmp(Species(:,3),'LN' ));...
    find(strcmp(Species(:,1),'mAPC')& strcmp(Species(:,3),'LN' ))];   
plot(simData.time/30, simData.data(:,i), 'LineWidth', linesize)
ylabel('Number of cells','Fontsize',fontsize); 
legend('Tumor APC','Tumor mAPC','TdLN APC','TdLN mAPC','Location','southeast')
set(gca, 'YScale', 'log'); ylim([1 1e9])   
xlabel('Time (month)','Fontsize',fontsize);
xlim([min(simData.time/30) max(simData.time/30)]);
set(gca,'TickLength',[0.03, 0.015])
set(gca,'XTick',[0, 6, 12])

subplot(3,2,3); hold on; box on; clear i; 
i =[find(strcmp(Species(:,1),'Cp')       & strcmp(Species(:,3),'Tum')); ...
    find(strcmp(Species(:,1),'D1_0')     & strcmp(Species(:,3),'Tum'))];   
plot(simData.time/30, simData.data(:,i)*1e6, 'LineWidth', linesize)
ylabel('Protein in Tumor (\muM)','Fontsize',fontsize);
legend('Self','Antigen','Location','southeast')
set(gca, 'YScale', 'log');  ylim([1e-3 1e0])   
xlabel('Time (month)','Fontsize',fontsize);
xlim([min(simData.time/30) max(simData.time/30)]);
set(gca,'TickLength',[0.03, 0.015])
set(gca,'XTick',[0, 6, 12])

subplot(3,2,4); hold on; box on; clear i; 
i =[find(strcmp(Species(:,1),'cpt_M')  & strcmp(Species(:,3),'mAPC_Surf'));...
    find(strcmp(Species(:,1),'p1_0_M') & strcmp(Species(:,3),'mAPC_Surf'))];
plot(simData.time/30, simData.data(:,i), 'LineWidth', linesize)
ylabel('Displayed molecules on mAPC','Fontsize',fontsize); 
legend('Self','Antigen','Location','southeast')
set(gca, 'YScale', 'log'); ylim([1e3 1e6])   
xlabel('Time (month)','Fontsize',fontsize);
xlim([min(simData.time/30) max(simData.time/30)]);
set(gca,'TickLength',[0.03, 0.015])
set(gca,'XTick',[0, 6, 12])

subplot(3,2,5); hold on; box on; clear i; 
i =[find(strcmp(Species(:,1),'Nivo')& strcmp(Species(:,3),'Cent')), ...
    find(strcmp(Species(:,1),'Nivo')& strcmp(Species(:,3),'Tum'))];
f =[model.Parameters(find(strcmp(Parameters(:,1),'f_vol_cent_Nivo' ))).Value, ...
    model.Parameters(find(strcmp(Parameters(:,1),'f_vol_tum'       ))).Value]; 
MW_Nivo = model.Parameters(find(strcmp(Parameters(:,1),'MW_Nivo' ))).Value;
plot(simData.time/30, simData.data(:,i(1))*MW_Nivo./f(1), 'LineWidth', linesize)
plot(simData.time/30, simData.data(:,i(2))*MW_Nivo./f(2), 'LineWidth', linesize)
legend('Central','Tumor','Location','northeast')
ylabel('Nivolumab (\mug/mL)','Fontsize',fontsize);
xlabel('Time (month)','Fontsize',fontsize);
xlim([min(simData.time/30) max(simData.time/30)]);
set(gca,'TickLength',[0.03, 0.015])
set(gca,'XTick',[0, 6, 12])

subplot(3,2,6); hold on; box on; clear i; 
i =[find(strcmp(Species(:,1),'C1_PDL1_Teff_PD1') & strcmp(Species(:,3),'Tum')); ...
    find(strcmp(Species(:,1),'C1_PDL2_Teff_PD1') & strcmp(Species(:,3),'Tum'))];   
plot(simData.time/30, simData.data(:,i), 'LineWidth', linesize)
ylabel('Number of molecules in synapse','Fontsize',fontsize);
legend('PDL1_{C}-PD1_{Teff}','PDL2_{C}-PD1_{Teff}','Location','northeast')
set(gca, 'YScale', 'log'); ylim([1e-0 1e5]);
xlabel('Time (month)','Fontsize',fontsize);
xlim([min(simData.time/30) max(simData.time/30)]);
set(gca,'TickLength',[0.03, 0.015])
set(gca,'XTick',[0, 6, 12])