
addpath DREAM5_NetworkInferenceChallenge_AlternativeDataFormats/net1/
addpath Network_predictions/Network_predictions/Community' 'integration

CCC=importdata('DREAM5_NetworkInference_Community_Network1.txt');

BBB=importdata('net1_expression_data_avg_t.tsv','\t');

%AAA=xlsread('DREAM5_NetworkInference_Ecoli_Saureus_Modules.xls');