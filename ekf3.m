%Load Files
clear
close
load dataset2.mat

%Variances values
T = 0.1;
%Q_k = [v_var, 0, 0;0, v_var,0;0, 0, om_var];
Q_k = diag([v_var;om_var])*10;

%Initialization of arrays
sizet = size(t);
iterations = 1000;%sizet(1,1);
F_kprev = [];
w_kprev = [];
x_kprev = zeros(iterations,3);
p_kprev = [];
x_kpost = zeros(iterations,3);
p_kpost = [];
r_thresh = 1;
variances = zeros(3,iterations);
cov_matrix = zeros(iterations,4);

for i = 1:iterations
    
    if i ==1
        
        % Prev Covariance matrix
        F_kprev = [1, 0 , -T*sin(th_true(1))*v(i); 0,1, T*cos(th_true(1))*v(i); 0,0,1];
        w_kprev = [cos(th_true(1)),0;sin(th_true(1)),0;0,1]*T;
        p_kprev = F_kprev * (diag([1,1,0.1])) * F_kprev' + w_kprev * Q_k * w_kprev';
   
        % Prev State
        vec1 = [x_true(1); y_true(1); th_true(1)] + T*[cos(th_true(1)),0;sin(th_true(1)),0;0,1] * [v(i);om(i)];
        x_kprev(i,1:3) = vec1';
        x_kprev(i,3) = wrapToPi(x_kprev(i,3));
        
        a = r(i,:);
        a = a(a~=0);
        
        a = a(:, ~(any(a > r_thresh,1)));
        
        if isempty(a) == 1  
           
            x_kpost(i,1:3) = x_kprev(i,1:3);
            p_kpost = p_kprev;
            variances(1:3,i) = diag(p_kpost);
            cov_matrix(i,1:2) = p_kpost(1,1:2);
            cov_matrix(i,3:4) = p_kpost(2,1:2);

            
        else
            
            %Kalman Gain & correction
            [K, G, y_kmeas, g] = kalmanGain2(x_kprev(i,1:3),p_kprev, r(i,:), b(i,:), r_thresh, l, d,r_var,b_var);
            %Posterior state
            vec2 = x_kprev(i,1:3)' + K*(y_kmeas - g);
            x_kpost(i,1:3) = vec2';
            x_kpost(i,3) = wrapToPi(x_kpost(i,3));
            %Posterior covariance
            p_kpost = (eye(3)-K*G)*p_kprev;
            variances(1:3,i) = diag(p_kpost);
            cov_matrix(i,1:2) = p_kpost(1,1:2);
            cov_matrix(i,3:4) = p_kpost(2,1:2);
            
        end
        
    else
        % Prev Covariance matrix
        F_kprev = [1, 0 , -T*sin(x_kpost(i-1,3))*v(i); 0,1, T*cos(x_kpost(i-1,3))*v(i); 0,0,1];
        w_kprev = [cos(x_kpost(i-1,3)),0;sin(x_kpost(i-1,3)),0;0,1]*T;
        p_kprev = F_kprev * p_kpost * F_kprev' + w_kprev * Q_k * w_kprev';
              
   
        % Prev State
        vec1 = [x_kpost(i-1,1); x_kpost(i-1,2); x_kpost(i-1,3)] + T*[cos(x_kpost(i-1,3)),0;sin(x_kpost(i-1,3)),0;0,1] * [v(i);om(i)];
        x_kprev(i,1:3) = vec1';
        x_kprev(i,3) = wrapToPi(x_kprev(i,3));
        
        a = r(i,:);
        a = a(a~=0);
        
        a = a(:, ~(any(a > r_thresh,1)));
        
        if isempty(a) == 1
           
            x_kpost(i,1:3) = x_kprev(i,1:3);
            p_kpost = p_kprev;
            variances(1:3,i) = diag(p_kpost);
            cov_matrix(i,1:2) = p_kpost(1,1:2);
            cov_matrix(i,3:4) = p_kpost(2,1:2);
            
        else
            
            %Kalman Gain & correction
            [K, G, y_kmeas, g] = kalmanGain2(x_kprev(i,1:3),p_kprev, r(i,:), b(i,:), r_thresh, l, d,r_var,b_var);
            %Posterior state
            vec2 = x_kprev(i,1:3)' + K*(y_kmeas - g);
            x_kpost(i,1:3) = vec2';
            x_kpost(i,3) = wrapToPi(x_kpost(i,3));
            %Posterior covariance
            p_kpost = (eye(3)-K*G)*p_kprev;
            variances(1:3,i) = diag(p_kpost);
            cov_matrix(i,1:2) = p_kpost(1,1:2);
            cov_matrix(i,3:4) = p_kpost(2,1:2);
            
        end
  
    end
end