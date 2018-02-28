function  [ModelInfo] = Kriging_train(S, Y)

[m, n] = size(S);
lob = repmat(1e-3,1,n);
upb = repmat(100,1,n);

%% Check design points
global ModelInfo;
[ModelInfo.theta, MinNegLnLikelihood] = ga(@likelihood, n,...
    [], [], [], [], lob, upb);
[NegLnLike, ModelInfo.R, ModelInfo.U] = likelihood(ModelInfo.theta);
end
