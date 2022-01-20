function newKalmanS = rtmyKalman(oldKalmanS,feats)

Y_t = feats';
Pk_apo = oldKalmanS.Pk_apo;
x_apo = oldKalmanS.x_apo;

numStateVar = oldKalmanS.numStateVar; 

A = oldKalmanS.A;
W = oldKalmanS.W;
C = oldKalmanS.C;
Q = oldKalmanS.Q;

newKalmanS = oldKalmanS;



%x_apr = 
% TIME UPDATE EQUATIONS
% STEP 1
x_apr =  A * x_apo;
% STEP 2
Pk_apr = A * Pk_apo * A' + W;


% MEASUREMENT UPDATE EQUATIONS
try
    % STEP 3 
    KK = Pk_apr * C' * pinv(C * Pk_apr * C' + Q); % Filter gain

    % STEP 4
    newKalmanS.x_apo = x_apr + KK * (Y_t - C * x_apr); % A posteriori Estimate
    
    % STEP 5
    newKalmanS.Pk_apo = (ones(numStateVar,numStateVar) - KK * C) * Pk_apr;
    
catch
    newKalmanS.Pk_apo = diag(1.*ones(1,numStateVar));
    newKalmanS.x_apo = x_apr;
end



end