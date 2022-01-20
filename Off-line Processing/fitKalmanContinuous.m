function myTrainedKalman = fitKalmanContinuous(X_t, Y_tr, nS)




%%
%%
%%
%% Traditional Kalman Filter approach
%%

%%%%***     Build the Kalman Filter Model 
%%%%***     X_t = A*X_(t-1) + W_t;
%%%%***     Y_t = C*X_t + Q_t;


% These are one-step shifted state matrices
X_1 = X_t(:,1:end-1);
X_2 = X_t(:,2:end);

% Calculate the coefficient matrices
% H = (Z_t*X_t')*inv(X_t*X_t'); % Slightly different nomenclature used elsewhere
% C = (Y_t * X_t') * (inv(X_t * X_t'));
C = (Y_tr * X_t') * (pinv(X_t * X_t'));

Q = (1/nS) .* ((Y_tr - (C * X_t)) * (Y_tr - (C * X_t))'); % Q = (Z_t-H*X_t)*(Z_t-H*X_t)' ./ (training_length);

% A = (X_2 * X_1') * inv(X_1 * X_1'); %  A = (X_2*X_1')*inv(X_1*X_1');
A = (X_2 * X_1') * pinv(X_1 * X_1'); %  A = (X_2*X_1')*inv(X_1*X_1');

W = (1/(nS-1)) .* (X_2 - A * X_1) * ((X_2 - A * X_1)');
% W = ((X_2-A*X_1)*(X_2-A*X_1)')./(training_length-1);

%% Let's run the prediction now
% Pk_apr = zeros(size(X_t,1),size(X_t,1));
% Pk_apo = ones(size(X_t,1),size(X_t,1));
Pk_apo = diag(1.*ones(1,size(X_t,1)));
% Pk_apo = diag(ones(1,size(X_t,1)));
x_apr = zeros(size(X_t,1),1);
x_apo = zeros(size(X_t,1),1);

myTrainedKalman.A = A;
myTrainedKalman.W = W;
myTrainedKalman.C = C;
myTrainedKalman.Q = Q;

myTrainedKalman.Pk_apo = Pk_apo;
myTrainedKalman.x_apo = x_apo;
myTrainedKalman.numStateVar = size(X_t,1);

