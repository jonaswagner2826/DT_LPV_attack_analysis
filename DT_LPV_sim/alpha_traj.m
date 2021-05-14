function [p] = alpha_traj(k)
    %alpha_traj Function returns the parameter as a function of k
    k_limit = 50;
    p1 = [0.50; 0.30; 0.20; 0.00];
    p2 = [0.35; 0.40; 0.10; 0.15];
    if (k > k_limit)
        p = p2;
    else
        p = p1;
    end
end

