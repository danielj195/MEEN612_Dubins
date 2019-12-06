function [ V_ordered ] = order_vertices( V, D )
%Inputs:  V-  nx2 matrix of vertices 
%         D-  order in which vertices are visited
%Outputs: V_ordered-  nx2 matrix of ordered verticies

V_temp = zeros(length(D),2);
for i = 1:length(D)
    V_temp(i,:) = V(D(i),:);
end
V_ordered = V_temp;
end

