predictionfile = '../INPUT/predictions/myteam/DREAM5_NetworkInference_myteam_Network1.txt';

M = load_dream_network(predictionfile);

s=uint16(M(:,1));
t=uint16(M(:,2));
G = digraph(s',t');

plot(G,'Layout','force')