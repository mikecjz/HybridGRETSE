fileID = fopen('/mnt/radnas1/Junzhou/SimulationProtocol.txt','r');
A = fscanf(fileID,'%c');
fclose(fileID);

A = A(300:end);

Numbers = split(A,' ');

Numbers = Numbers(~(cellfun('isempty',Numbers)));

RF = Numbers(2:2:end);
RF = str2double(RF);

RF_degree = RF * 74.6/257;

RF_degree = RF_degree(RF_degree~=0);

RF_degree = unique(RF_degree,'stable');