function  [K, G, y_kmeas, g] = kalmanGainCRLB(x_prev,x_true,p_prev, r_i, b_i, r_tresh, positions, dist,r_vari,b_vari)
%Compute the Kalman Gain K, G jacobian, and y_kmeas measurements available
%x_prev(1,3) = wrapToPi(x_prev(1,3));

%Eliminate zero range values
positions = positions';
colsWithZeros = any(r_i==0,1);
r_i = r_i(:, ~colsWithZeros);
b_i = b_i(:, ~colsWithZeros);
positions = positions(:, ~colsWithZeros);

%Eliminate range values higher than r_tresh
colsWithTresh = any(r_i > r_tresh,1);
r_i = r_i(:, ~colsWithTresh);
b_i = b_i(:, ~colsWithTresh);
positions = positions(:, ~colsWithTresh);

%Obtain size of measurements
sizeTotal = size(r_i);
sizeTotal = sizeTotal(1,2);

%Compute y_kmeas: stacked vectors of range and bearing
y_kmeas(1:sizeTotal,1) = r_i';
y_kmeas(sizeTotal+1:sizeTotal*2,1) = b_i';

%Computation of G
G = zeros(sizeTotal*2,3);
g = zeros(sizeTotal*2,1);
%Expressions for observation model
syms xl xk tetak d yl yk nk
R_K = sqrt((xl-xk-(d*cos(tetak)))^2 + (yl-yk-(d*sin(tetak)))^2) + nk;
fi_k = atan2(yl-yk-(d*sin(tetak)),xl-xk-(d*cos(tetak))) - tetak + nk;

G1= jacobian(R_K,[xk,yk,tetak]);
G2= jacobian(fi_k,[xk,yk,tetak]);



for i = 1:sizeTotal
    
    G(i,1:3) = vpa(subs(G1,[xk, yk, tetak, xl, yl, d, nk],[x_true(1,1), x_true(1,2), x_true(1,3), positions(1,i), positions(2,i), dist,0]));
    
end


for j = 1:sizeTotal
    
    G(sizeTotal+j,1:3) = vpa(subs(G2,[xk, yk, tetak, xl, yl, d, nk],[x_true(1,1), x_true(1,2), x_true(1,3), positions(1,j), positions(2,j), dist,0]));

end

%Then calculation of nk_prime, considering derivative is = 1
r_k = eye(sizeTotal*2);
r_k(1:sizeTotal,:) = r_k(1:sizeTotal,:)*r_vari;
r_k(sizeTotal+1:sizeTotal*2,:) = r_k(sizeTotal+1:sizeTotal*2,:)*b_vari;

%Finally, compute Kalman Gain K

K = (p_prev * G')/(G * p_prev * G' + r_k);

%And g values
for i = 1:sizeTotal
    
    g(i,1) = vpa(subs(R_K,[xk, yk, tetak, xl, yl, d, nk],[x_prev(1,1), x_prev(1,2), x_prev(1,3), positions(1,i), positions(2,i), dist,0]));
    
end


for i = 1:sizeTotal
    
    g(sizeTotal+i,1) = vpa(subs(fi_k,[xk, yk, tetak, xl, yl, d, nk],[x_prev(1,1), x_prev(1,2), x_prev(1,3), positions(1,i), positions(2,i), dist,0]));
    g(sizeTotal+i,1) = wrapToPi(g(sizeTotal+i,1));
    
end


end

