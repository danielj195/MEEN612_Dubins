function [ Cost ] = path_cost_matrix( V_list, angles, r )
%Inputs: V_list - list of vertices to be visited obtained by doubling tree
%               heuristic (nx2) matrix; n vertices each with x,y coordinate 
%       Angles - List of heading angles to be considered at each vertex
%               (1xa) vector
%       r-      minimum turning radius
%Outputs: Adj - adjacency matrix for all possible dubins paths (n*a x n*a)
%matrix
A_len = size(V_list,1)*length(angles);
A = zeros(A_len, A_len);

for i = 1:(size(V_list,1)-1)
    block = zeros(length(angles),length(angles));
    for j = 1:length(angles)  %initial heading angle
        for k = 1:length(angles)  %final heading angle
            p1 = [V_list(i,:), angles(j)];  %initial config
            p2 = [V_list(i+1,:), angles(k)];  %final config
            param = dubins_core(p1,p2,r);
            block(j,k) = dubins_length(param);
        end
    end
    %Insert block into adjacency matrix
    r_start = i*length(angles) - length(angles) + 1;
    r_end = i*length(angles);
    c_start = i*length(angles) + 1;
    c_end = i*length(angles) + length(angles);
    A(r_start:r_end, c_start:c_end) = block;
end

% %Connect final vertex to initial vertex
% block = zeros(length(angles),length(angles));
% for j = 1:length(angles)  %initial heading angle
%     for k = 1:length(angles)  %final heading angle
%         p1 = [V_list(size(V_list,1),:), angles(j)];  %initial config
%         p2 = [V_list(1,:), angles(k)];  %final config
%         param = dubins_core(p1,p2,r);
%         block(j,k) = dubins_length(param);
%     end
% end
% 
% %Insert block into adjacency matrix
% r_start = size(V_list,1)*length(angles) - length(angles) + 1;
% r_end = size(V_list,1)*length(angles);
% c_start = 1;
% c_end = length(angles);
% A(r_start:r_end, c_start:c_end) = block;

Cost = A;
end

