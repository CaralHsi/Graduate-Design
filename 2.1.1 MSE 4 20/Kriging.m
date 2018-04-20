function Kriging( sample_point, sample_value)
% input sample points and sample value and x to be predicted
% output the value at x
%% clear unused data
clear global;

%% load data
[sample_point,ia,ic] = unique(sample_point,'rows');
sample_value = sample_value(ia);
S = sample_point;
Y = sample_value;
global ModelInfo
ModelInfo.S_oringinal = S;
ModelInfo.Y_oringinal = Y;
m = size(S,1);  % number of design sites and their dimension

%% Normalize data
mS = mean(S);   sS = std(S);
mY = mean(Y);   sY = std(Y);
% 02.08.27: Check for 'missing dimension'
% %WHY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
j = find(sS == 0);
if  ~isempty(j),  sS(j) = 1; end
j = find(sY == 0);
if  ~isempty(j),  sY(j) = 1; end
S = (S - repmat(mS,m,1)) ./ repmat(sS,m,1);
Y = (Y - repmat(mY,m,1)) ./ repmat(sY,m,1);



%% set global variable
global ModelInfo
ModelInfo.S = S;
ModelInfo.Y = Y;
global Normalization
Normalization.mS = mS;
Normalization.sS = sS;
Normalization.mY = mY;
Normalization.sY = sY;
Kriging_train(S, Y);

end



