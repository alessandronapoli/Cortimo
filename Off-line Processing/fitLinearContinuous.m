function myLinearMod = fitLinearContinuous(X_t, Y_tr)




%%
%%
%% Linear Model
%%

linearmodX = fitlm(Y_tr',X_t(1,:)','linear');
linearmodY = fitlm(Y_tr',X_t(2,:)','linear');


myLinearMod.X = linearmodX;
myLinearMod.Y = linearmodY;


