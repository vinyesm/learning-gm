clear all; 
load ..\..\data\data.mat


%% -------------------------------- Greedy-Feige for finding multiple DkS -------------------------------
  
dim = size(W,1);

cardinality_sequence = [30 30 30 30 30 30];

idx = 1:dim;
W_cur = W;
dense_subgraph_idx_mul = [];
center_city_idx_mul = [];


% Output Log
fprintf('Greedy-Feige %10s %10s %10s     %10s \n', ' Subgraph ID ', ' Cardinality ', ' Density ', ' CPU Time ');

% Finding multiple DkS
for (i = 1:length(cardinality_sequence))
    
    k = cardinality_sequence(i);

    tic;
    u = Greedy_DkS_Feige(W_cur, k); 
    t_cpu = toc;
    Density =  u'*W_cur*u / k;
    
    fprintf('Greedy-Feige %10d   %10d     %10f     %10f  \n', i, k, Density, t_cpu);
        
    idx_1 = idx(find(u~=0));
    idx = idx(find(u==0));
    W_cur = W(idx, idx);
    
    city_name(idx_1)';
    
    deg = sum(W(idx_1, idx_1),2);
    [val , id] = max(deg);
    center_city_idx_mul= [center_city_idx_mul, idx_1(id)];
    
    dense_subgraph_idx_mul = [dense_subgraph_idx_mul, {idx_1}];

end

% Save results
Greedy_Feige_Multiple_Result = [];
Greedy_Feige_Multiple_Result.cardinality = cardinality_sequence;
Greedy_Feige_Multiple_Result.dense_subgraph_idx = dense_subgraph_idx_mul;
Greedy_Feige_Multiple_Result.center_city_idx = center_city_idx_mul;

save ..\..\result\Greedy_Feige_Multiple_Result.mat 'Greedy_Feige_Multiple_Result';

