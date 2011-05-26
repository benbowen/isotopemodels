%a script to set values for all required parameters
function params = set_nuclei_composition(params)
% clear params;
% 
% params.charge = 1;
% params.resolution = 0.1;
% params.delta_m = 5e-6;
% params.delta_i = 0.1;
% params.numnoisepeaks = 20;
% params.noisepeakintensity = 1e3;
% params.I_max = 1e8;
% params.I_min = 1e2;
% params.model = 0;

%%
load Fractionation_El_Mass_Abund.mat

% padd it with 20 Hn# for the global fitting analysis
s=size(mass,1);
for i = 1:20
    ele=cat(2,ele,['Hn',num2str(i)]);
    mass = cat(1,mass,mass(s,:));
    abund = cat(1,abund,abund(s,:));
end


params.nucleiinfo.ele = ele;
params.nucleiinfo.mass = mass;
params.nucleiinfo.abund = abund;

clear ele mass abund;