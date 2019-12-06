function [ A ] = complete_graph( V )
%Inputs:  V- nx2 matrix where row n is planar coordinate of vertex n
%Outputs: A- Adjacency matrix where element (i,j) is euclidean distance
%from vertex i to vertex j;  Matrix rows/columns are in same order as the
%vertices in matrix V
A = zeros(length(V),length(V));
for i = 1:length(V)  
    for j = 1:length(V)   
        X = [V(i,:);V(j,:)];
        dist = pdist(X,'euclidean');
        A(i,j) = dist;
    end
end
end

