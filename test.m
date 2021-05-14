A(:,:,1) = [-0.80, 0.25; 0.25,-0.30]; B(:,:,1) = [ 1.90; 0.00];
A(:,:,2) = [ 0.30, 0.70; 0.70, 0.00]; B(:,:,2) = [-1.00; 1.50];
A(:,:,3) = [-0.30, 0.65; 0.55, 0.10]; B(:,:,3) = [ 0.30;-2.00];
A(:,:,4) = [ 0.55,-0.20;-0.40,-0.30]; B(:,:,4) = [-0.60; 0.00];
C = [1 0];


Alpha = sym('alpha',[4,1]);

n = size(A,1);
m = size(A,3);
p = size(B,2);
q = size(C,1);

A_bar = zeros(n);
B_bar = zeros(n,p);
for i =1:m
    A_bar = A_bar + Alpha(i) * A(:,:,i);
    B_bar = B_bar + Alpha(i) * B(:,:,i);
end

syms q
a_coeff = det(q*eye(n) - A_bar);
charpoly(a_coeff,q);

% - q^2
% q*((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1)
% (11*alpha1*alpha2)/25 + (29*alpha1*alpha3)/100 - (9*alpha1*alpha4)/40 + (81*alpha2*alpha3)/100 - (33*alpha2*alpha4)/100 - (103*alpha3*alpha4)/200 - (71*alpha1^2)/400 + (49*alpha2^2)/100 + (31*alpha3^2)/80 + (49*alpha4^2)/200

b_coeff = C * adjoint(q*eye(n) - A_bar) * B_bar;
charpoly(b_coeff,q)

% 0
% q*(alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1)
% (3*alpha1*alpha3)/5 - (3*alpha1*alpha2)/40 - (39*alpha1*alpha4)/100 + (13*alpha2*alpha3)/40 + (3*alpha2*alpha4)/5 - (11*alpha3*alpha4)/20 - (57*alpha1^2)/100 - (21*alpha2^2)/20 + (133*alpha3^2)/100 + (9*alpha4^2)/50


attempt1 = [
    ((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1);
    (11*alpha1*alpha2)/25 + (29*alpha1*alpha3)/100 - (9*alpha1*alpha4)/40 + (81*alpha2*alpha3)/100 - (33*alpha2*alpha4)/100 - (103*alpha3*alpha4)/200 - (71*alpha1^2)/400 + (49*alpha2^2)/100 + (31*alpha3^2)/80 + (49*alpha4^2)/200;
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1);
    (3*alpha1*alpha3)/5 - (3*alpha1*alpha2)/40 - (39*alpha1*alpha4)/100 + (13*alpha2*alpha3)/40 + (3*alpha2*alpha4)/5 - (11*alpha3*alpha4)/20 - (57*alpha1^2)/100 - (21*alpha2^2)/20 + (133*alpha3^2)/100 + (9*alpha4^2)/50;
    ];
attempt2 = [
    (11*alpha1*alpha2)/25 + (29*alpha1*alpha3)/100 - (9*alpha1*alpha4)/40 + (81*alpha2*alpha3)/100 - (33*alpha2*alpha4)/100 - (103*alpha3*alpha4)/200 - (71*alpha1^2)/400 + (49*alpha2^2)/100 + (31*alpha3^2)/80 + (49*alpha4^2)/200;
    ((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1);
    (3*alpha1*alpha3)/5 - (3*alpha1*alpha2)/40 - (39*alpha1*alpha4)/100 + (13*alpha2*alpha3)/40 + (3*alpha2*alpha4)/5 - (11*alpha3*alpha4)/20 - (57*alpha1^2)/100 - (21*alpha2^2)/20 + (133*alpha3^2)/100 + (9*alpha4^2)/50;
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1)
    ];
attempt3 = [
    ((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1);
    -1;
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1);
    0
    ];
attempt4 = [
    -1;
    ((3*alpha2)/10 - (11*alpha1)/10 - alpha3/5 + alpha4/4 + 1);
    0;
    (alpha2 - (19*alpha1)/10 - (3*alpha3)/10 + (3*alpha4)/5 + 1)
    ];