function [x_hat, alpha_hat, nu_hat] ...
                                        = est_DE(x_hat, alpha_hat, y, u,...
                                           P, phi, nu_hat,...
                                           gamma, A, B, C)
    %est_JSPE Function performs a single iteration of the joint state and
    %parameter estimator
    
    arguments
        x_hat       %x_hat_{k-1}        (n,1)
        alpha_hat   %\alpha_hat_{k-1}   (m,1)
        y           %y_k                (1,1)
        u           %u_k                (1,1)
        P           %P_k                (2n,2n)
        phi         %phi_k              (2n,1)
        nu_hat      %nu_hat_k           (1,2n)
        gamma = 0.9
        A = string(-1)
        B = string(-1)
        C = [1, 0]
%         nu_alpha = -1
    end
    
    %% System Definition
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
    
    % System Dimenstions
    n = size(A,1);
    m = size(A,3);
    p = size(B,2);
    q = size(C,1);
    
    %% State Estimation
    % CVX Feasability Problem on LMI
    disp('DE State Estimation Started')
    tol = 1e-6;
    cvx_begin sdp quiet
        variable P_cvx(n,n,m) symmetric
        variable F(n,q,m)
        variable G(n,n,m)
        variable zeta
        subject to
            for i = 1:m
                for j = 1:m
                    [G(:,:,i) + G(:,:,i)' - P_cvx(:,:,j), zeros(n), G(:,:,i)*A(:,:,i) + F(:,:,i)*C, G(:,:,i);
                     zeros(n), eye(n), eye(n), zeros(n);
                     A(:,:,i)'*G(:,:,i)' + C'*F(:,:,i)', eye(n), P_cvx(:,:,i), zeros(n);
                     G(:,:,i)', zeros(n), zeros(n), zeta*eye(n)] >= tol*eye(4*n);
                end
            end
            zeta >= 1;
    cvx_end

    % Gain Result Calc
    for i = 1:m
        L(:,:,i) = inv(G(:,:,i))*F(:,:,i);
    end
    
    % State Estimate Update
    x_hat_old = x_hat; % x_hat_{k-1}
    x_hat = 0;
    for i = 1:m
        x_hat = x_hat + alpha_hat(i)*(A(:,:,i)*x_hat_old + B(:,i)*u)...
            + L(:,:,i)*(C*x_hat_old - y);
    end
    %Note: x_hat = x_hat k
    
    %% Parameter Estimation
   
%     % Phi Matrix Definitions
%     if string(nu_alpha) == string(-1)
%         Alpha_sym = sym('alpha',[m,1]);
%         A_bar = zeros(n);
%         B_bar = zeros(n,p);
%         for i = 1:m
%             A_bar = A_bar + Alpha_sym(i) * A(:,:,i);
%             B_bar = B_bar + Alpha_sym(i) * B(:,:,i);
%         end
%         a_coeff = charpoly(A_bar);
%         temp = sym('temp');
%         b_coeff = charpoly(C*adjoint(temp*eye(n) - A_bar)*B_bar);
%         if size(a_coeff,2) <= n
%             a_coeff = [zeros(1,n+1-size(a_coeff,2)),a_coeff];
%         end
%         if size(b_coeff,2) <= n
%             b_coeff = [zeros(1,n+1-size(b_coeff,2)),b_coeff];
%         end
%         nu_alpha = matlabFunction([fliplr(a_coeff(1:n)), fliplr(b_coeff(1:n))]');
%     end
    

%     % Est. Opt. Problem
%     disp('DE Parameter Estimation Started')
%     % attempting yalmip
%     yalmip('clear')
%     Alpha = sdpvar(m,1);
% %     nu_alpha = sdpvar(2*n,1);
%     
%     Constraints = [sum(Alpha) == 1];
% %     Constraints = [Constraints, nu_alpha == nuVector(Alpha)];
%     for i = 1:m
%         Constraints  = [Constraints, Alpha(i) >= 0];
%     end
%     
%     Objective = (nu_hat - nu_alpha)' * inv(P) * (nu_hat - nu_alpha);
%     
%     options = sdpsettings();%'debug',1,'verbose',1)%,'solver','quadprog','quadprog.maxiter',100);
%     
%     sol = optimize(Constraints,Objective,options)
%     
%     
% %     cvx_begin quiet
% %         variable alpha_cvx(m,1)
% %         nu_alpha = nu_alpha(alpha_cvx);
% %         
% %         % Minimize This
% %         minimize((nu_hat - nu_alpha)' * inv(P) * (nu_hat - nu_alpha))
% %         subject to
% %             sum = 0;
% %             for i = 1:m
% %                 alpha_cvx(i) >= 0;
% %                 sum = sum + alpha_cvx(i);
% %             end
% %             sum == 1;
% %     cvx_end
%     if sol.problem == 0
%          % Extract and display value
%          alpha_hat = value(Alpha) %alpha_cvx;
%          nu_hat = value(nu_alpha)
%     else
%          disp('Hmm, something went wrong!');
%          sol.info
%          yalmiperror(sol.problem)
%          error('issue here.....')
%     end
    
    
    
%     % k step
%     P_old = P;
%     phi_old = phi;
%     nu_hat_old = nu_hat;
% 
% 
% %     % k+1 step
% %     phi(2:n) = phi_old(1:n-1);
% %     phi(1) = y;
% %     phi(n+2:2*n) = phi_old(n+1:2*n-1);
% %     phi(n+1) = u;
%     
%     % Recursive Update
%     P = (1/gamma)*P_old - (1/gamma)*(phi'...
%         * inv(gamma + phi * P_old * phi') * phi * P_old);
%     nu_hat = nu_hat_old + P * phi' * (y - phi*nu_hat_old);    
    
    

    
    
end

