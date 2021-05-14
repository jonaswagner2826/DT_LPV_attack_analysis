function nu = nuVector(Alpha,A,B,C)
    %UNTITLED8 Summary of this function goes here
    %   Detailed explanation goes here
    arguments
        Alpha
        A = 1;
        B = 1;
        C = 1;
    end
    
%     nu = [
%         
%         11*alpha(1)/10 - 3*alpha(2)/10 + alpha(3)/5 - alpha(4)/4;
% %         1;
%         1;
%         0];

%     [alpha1, alpha2, alpha3, alpha4] = Alpha;
    alpha1 = Alpha(1);
    alpha2 = Alpha(2);
    alpha3 = Alpha(3);
    alpha4 = Alpha(4);
%     nu = [
% (71*(alpha1)^2)/400 - (49*(alpha2)^2)/100 - (31*(alpha3)^2)/80 - (49*(alpha4)^2)/200 - (11*(alpha1)*(alpha2))/25 - (29*(alpha1)*(alpha3))/100 + (9*(alpha1)*(alpha4))/40 - (81*(alpha2)*(alpha3))/100 + (33*(alpha2)*(alpha4))/100 + (103*(alpha3)*(alpha4))/200;
% (11*(alpha1))/10 - (3*(alpha2))/10 + (alpha3)/5 - (alpha4)/4;
% (133*(alpha3)^2)/100 - (21*(alpha2)^2)/20 - (57*(alpha1)^2)/100 + (9*(alpha4)^2)/50 - (3*(alpha1)*(alpha2))/40 + (3*(alpha1)*(alpha3))/5 - (39*(alpha1)*(alpha4))/100 + (13*(alpha2)*(alpha3))/40 + (3*(alpha2)*(alpha4))/5 - (11*(alpha3)*(alpha4))/20 - (19*(alpha1)*(temp))/10 + (alpha2)*(temp) - (3*(alpha3)*(temp))/10 + (3*(alpha4)*(temp))/5;
%    1];
    
    
%     
%     n = size(A,1);
%     m = size(A,3);
%     p = size(B,2);
%     q = size(C,1);
%     
%     Alpha_sym = sym('alpha',[m,1]);
%     A_bar = zeros(n);
%     B_bar = zeros(n,p);
%     for i = 1:m
%         A_bar = A_bar + Alpha_sym(i) * A(:,:,i);
%         B_bar = B_bar + Alpha_sym(i) * B(:,:,i);
%     end
%     a_coeff = charpoly(A_bar);
%     temp = sym('temp');
%     b_coeff = charpoly(C*adjoint(temp*eye(n) - A_bar)*B_bar);
%     if size(a_coeff,2) <= n
%         a_coeff = [zeros(1,n+1-size(a_coeff,2)),a_coeff];
%     end
%     if size(b_coeff,2) <= n
%         b_coeff = [zeros(1,n+1-size(b_coeff,2)),b_coeff];
%     end
%     
% %     nu_alpha = matlabFunction([fliplr(a_coeff(2:n+1)), fliplr(b_coeff(2:n+1))]');
%     
%     nu = [fliplr(a_coeff(2:n+1)), fliplr(b_coeff(2:n+1))]';



attempt1 = [
    -((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1),
    -(11*alpha1*alpha2)/25 + (29*alpha1*alpha3)/100 - (9*alpha1*alpha4)/40 + (81*alpha2*alpha3)/100 - (33*alpha2*alpha4)/100 - (103*alpha3*alpha4)/200 - (71*alpha1^2)/400 + (49*alpha2^2)/100 + (31*alpha3^2)/80 + (49*alpha4^2)/200,
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1),
    (3*alpha1*alpha3)/5 - (3*alpha1*alpha2)/40 - (39*alpha1*alpha4)/100 + (13*alpha2*alpha3)/40 + (3*alpha2*alpha4)/5 - (11*alpha3*alpha4)/20 - (57*alpha1^2)/100 - (21*alpha2^2)/20 + (133*alpha3^2)/100 + (9*alpha4^2)/50
    ];
attempt2 = [
    -(11*alpha1*alpha2)/25 + (29*alpha1*alpha3)/100 - (9*alpha1*alpha4)/40 + (81*alpha2*alpha3)/100 - (33*alpha2*alpha4)/100 - (103*alpha3*alpha4)/200 - (71*alpha1^2)/400 + (49*alpha2^2)/100 + (31*alpha3^2)/80 + (49*alpha4^2)/200,
    -((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1),
    (3*alpha1*alpha3)/5 - (3*alpha1*alpha2)/40 - (39*alpha1*alpha4)/100 + (13*alpha2*alpha3)/40 + (3*alpha2*alpha4)/5 - (11*alpha3*alpha4)/20 - (57*alpha1^2)/100 - (21*alpha2^2)/20 + (133*alpha3^2)/100 + (9*alpha4^2)/50,
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1)
    ];
attempt3 = [
    -((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1),
    1,
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1),
    0
    ];
attempt4 = [
    1,
    -((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1),
    0,
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1)
    ];

nu = attempt4;

end

