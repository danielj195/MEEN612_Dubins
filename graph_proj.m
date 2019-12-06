%List of (x,y) coordinates
V = [1.5 3.5;
     4.2 2.7;
     0.2 3.8;
     2.5 1.3;
     3.7 4.7;
     2.1 0.2;
     0.3 0.8;
     1.7 4.6;
     4.8 1.2;
     2.4 1.0];
 
% V = [0 0;
%      1 0;
%      1 1;
%      0 1];
r = 0.2;  %minimum turning radius
stepsize = 0.05;  %Step size for plotting


A = complete_graph(V);
MST = prims_MST(A);
D = dtree_heuristic(A,MST);  %Doubling Tree Heuristic

D_short = D;
D_short(end) = [];

V_ordered = order_vertices(V,D); 


angles = [0 pi/2 pi 3*pi/2];
C = path_cost_matrix(V_ordered,angles,r); %cost matrix
Adj = double(C >0);                         %adjacency matrix

%Run Dijkstra's on first vertex for all heading angles
%Find path with lowest cost
cost_min = inf;
opt_nodes = [0;0];

for i=1:length(angles)
    start_node = i;
    end_node = size(V_ordered,1)*length(angles) - length(angles) + i;
    [cost,path_angles] = dijkstra(Adj,C,start_node,end_node);
    if cost < cost_min
        cost_min = cost;
%         opt_nodes = [start_node; end_node];
        opt_path = path_angles;
    end
end

cost_min
% p1 = [V_ordered(1,:),angles(4)];
% p2 = [V_ordered(2,:),angles(3)];
% path1 = dubins_curve(p1,p2,r,0.05,0);
% success = construct_path(path,V_ordered,angles,r)
% [cost,path] = dijkstra(Adj,C,2,18)
% [cost,path] = dijkstra(Adj,C,3,19)
% [cost,path] = dijkstra(Adj,C,4,20)


opt_path = mod(opt_path,length(angles));
opt_path(opt_path == 0) = length(angles); %replace zeros with last angle option



figure('name','Dubins curve');

for i = 1:(size(V_ordered,1)-1)
    p1 = [V_ordered(i,:),angles(opt_path(i))];
    p2 = [V_ordered(i+1,:),angles(opt_path(i+1))];
    param = dubins_core(p1, p2, r);
    path = dubins_path_sample_many(param, stepsize);

%     disp('dubins calculation time'); toc;
    % plotting
    tic;    % most of the time is spent on plotting
%     figure('name','Dubins curve');
    plot(path(:,1), path(:,2)); axis equal; hold on
    scatter(p1(1), p1(2), 45, '*','r','LineWidth',1); hold on;
    scatter(p2(1), p2(2), 45, 'square','b','LineWidth',1); hold on;
%     text(p1(1), p1(2),'start','HorizontalAlignment','center');
%     text(p2(1), p2(2),'end','VerticalAlignment','top');
%     disp('plot drawing time'); toc;
    hold on
end


function path = dubins_path_sample_many( param, stepsize)
    if param.flag < 0
        path = 0;
        return
    end
    length = dubins_length(param);
    path = -1 * ones(floor(length/stepsize), 3);
    x = 0;
    i = 1;
    while x <= length
        path(i, :) = dubins_path_sample( param, x );
        x = x + stepsize;
        i = i + 1;
    end
    return
end

function end_pt = dubins_path_sample(param, t)
    if( t < 0 || t >= dubins_length(param) || param.flag < 0)
        end_pt = -1;
        return;
    end

    % tprime is the normalised variant of the parameter t
    tprime = t / param.r;

    % In order to take rho != 1 into account this function needs to be more complex
    % than it would be otherwise. The transformation is done in five stages.
    %
    % 1. translate the components of the initial configuration to the origin
    % 2. generate the target configuration
    % 3. transform the target configuration
    %      scale the target configuration
    %      translate the target configration back to the original starting point
    %      normalise the target configurations angular component

    % The translated initial configuration
    p_init = [0, 0, param.p_init(3) ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The three segment types a path can be made up of
    L_SEG = 1;
    S_SEG = 2;
    R_SEG = 3;

    % The segment types for each of the Path types
    DIRDATA = [ L_SEG, S_SEG, L_SEG ;...
                L_SEG, S_SEG, R_SEG ;...
                R_SEG, S_SEG, L_SEG ;...
                R_SEG, S_SEG, R_SEG ;...
                R_SEG, L_SEG, R_SEG ;...
                L_SEG, R_SEG, L_SEG ]; 
    %%%%%%%%%%%%%%%%%%%%%%%%% END DEFINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate the target configuration
    types = DIRDATA(param.type, :);
    param1 = param.seg_param(1);
    param2 = param.seg_param(2);
    mid_pt1 = dubins_segment( param1, p_init, types(1) );
    mid_pt2 = dubins_segment( param2, mid_pt1,  types(2) );
    
    % Actual calculation of the position of tprime within the curve
    if( tprime < param1 ) 
        end_pt = dubins_segment( tprime, p_init,  types(1) );
    elseif( tprime < (param1+param2) ) 
        end_pt = dubins_segment( tprime-param1, mid_pt1,  types(2) );
    else 
        end_pt = dubins_segment( tprime-param1-param2, mid_pt2,  types(3) );
    end

    % scale the target configuration, translate back to the original starting point
    end_pt(1) = end_pt(1) * param.r + param.p_init(1);
    end_pt(2) = end_pt(2) * param.r + param.p_init(2);
    end_pt(3) = mod(end_pt(3), 2*pi);
    return;
end

function seg_end = dubins_segment(seg_param, seg_init, seg_type)
    L_SEG = 1;
    S_SEG = 2;
    R_SEG = 3;
    if( seg_type == L_SEG ) 
        seg_end(1) = seg_init(1) + sin(seg_init(3)+seg_param) - sin(seg_init(3));
        seg_end(2) = seg_init(2) - cos(seg_init(3)+seg_param) + cos(seg_init(3));
        seg_end(3) = seg_init(3) + seg_param;
    elseif( seg_type == R_SEG )
        seg_end(1) = seg_init(1) - sin(seg_init(3)-seg_param) + sin(seg_init(3));
        seg_end(2) = seg_init(2) + cos(seg_init(3)-seg_param) - cos(seg_init(3));
        seg_end(3) = seg_init(3) - seg_param;
    elseif( seg_type == S_SEG ) 
        seg_end(1) = seg_init(1) + cos(seg_init(3)) * seg_param;
        seg_end(2) = seg_init(2) + sin(seg_init(3)) * seg_param;
        seg_end(3) = seg_init(3);
    end
end