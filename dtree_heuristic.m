function [ D ] = dtree_heuristic(G, MST)
%Inputs:  G- Adjacency matrix (doubling matrix should be same size)
%         MST- minimum spanning tree
%Outputs:  D- Doubling Tree heuristic; 1xn array showing order that
%vertices should be visited

dG = zeros(length(G),length(G));
for i = 1:size(MST,1)
    dG(MST(i,1),MST(i,2)) = 1;
    dG(MST(i,2),MST(i,1)) = 1;
end
dG
vertex = 1;
tour_d = vertex;
a=1;
for i = 1:2*size(MST,1)
    a;
    v_list = find(dG(vertex,:));
    non_rep = setdiff(v_list,tour_d);
    if ~isempty(non_rep)
        next_vertex = non_rep(end);
    else
        next_vertex = v_list(end);
    end
%     next_vertex = v_list(end)
    dG(vertex,next_vertex) = 0;
    vertex = next_vertex;
    tour_d = [tour_d; vertex];
end
% tour_d

for i= 2:length(tour_d)
    if (ismember(tour_d(i),tour_d(1:(i-1))))
        tour_d(i) = 0;
    end
end
Dp = tour_d(tour_d~=0);
D = [Dp;tour_d(1)];
end

