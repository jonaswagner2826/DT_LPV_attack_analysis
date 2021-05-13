function [x_hat, alpha_hat, P] = est_EKF(x_hat, alpha_hat, P, y, u,...
                                         A, B, C, Q, R)
    %est_EKF Function performs a single iteration of the parameter and
    %state estimation using an EKF
    
    arguments
        x_hat
        alpha_hat
        P
        y
        u
        A = string(-1)
        B = string(-1)
        C = [1, 0]
        Q = string(-1)
        R = string(-1)
    end
    
    % System Matrices
    if string(A) == string(-1)
        clear A
        A(:,:,1) = [-0.80, 0.25; 0.25,-0.30];
        A(:,:,2) = [ 0.30, 0.70; 0.70, 0.00];
        A(:,:,3) = [-0.30, 0.65; 0.55, 0.10];
        A(:,:,4) = [ 0.55,-0.20;-0.40,-0.30];
    end
    if string(B) == string(-1)
        clear B
        B(:,:,1) = [ 1.90; 0.00];
        B(:,:,2) = [-1.00; 1.50];
        B(:,:,3) = [ 0.30;-2.00];
        B(:,:,4) = [-0.60; 0.00];
    end
    
    % System Dimenstion Matrices
    n = size(A,1);
    m = size(A,3);
    p = size(B,1);
    q = size(C,1);
    
    
    % EKF Covariance Matrices
    if string(Q) == string(-1)
        clear Q
        Q = diag([zeros(1,n), 100*ones(1,m)]);
    end
    if string(R) == string(-1)
        clear R
        R = 0.01;
    end
    
    % A_hat (k-1) calc
    A_hat = diag([zeros(1,n), ones(1,m)]);
    for i = 1:m
        A_hat(1:n,1:n) = A_hat(1:n,1:n) + alpha_hat(i) * A(:,:,i);
        A_hat(1:n,n+i) = A(:,:,i) * x_hat + B(:,:,i) * u;
    end
    
    
    % Preditiction Step
    x_hat_pre = 0;
    for i = 1:m
        x_hat_pre = x_hat_pre + alpha_hat(i) * (A(:,:,i) * x_hat + B(:,:,i) * u);
    end
    alpha_hat_pre = alpha_hat;
    P_pre = A_hat * P * A_hat' + Q;
    
    % L (k) calc
    C_hat = [C, zeros(1,m)];
    L = P_pre * C_hat' * inv(R + C_hat * P_pre * C_hat');
    
    % Update Step
    x_alpha_pre = [x_hat_pre; alpha_hat_pre];
    x_alpha_post = x_alpha_pre + L * (y - C * x_hat_pre);
    P_post = (eye(n+m) - L * C_hat) * P_pre;
    
    % Estimates
    x_hat = x_alpha_post(1:n);
    alpha_hat = x_alpha_post(n+1:n+m);
    P = P_post;
end

