function [forwBackVel, leftRightVel, rotVel, finish] = wavefront_main(...
    pts, contacts, position, orientation, varargin)

% Check varargin
if length(varargin) == 2
    goal_pos = varargin{1}; % Goal position
    env_map = varargin{2};  % Environment
    debug = true;           % Debug plot, it is true by default and can be specified as vararg{3}
elseif length(varargin) == 3
    goal_pos = varargin{1};
    env_map = varargin{2};
    debug = varargin{3};
else
    error("Wrong number of additional arguments provided: %d\n", length(varargin))
end

% get the current position and orientation
pos = position(1:2);
phi = orientation(3);

% Define max velocities and regulator gains
u_max = 10;
w_max = 10;
Kp_linear = 20;
Kp_angular = 10;

% Tolerance to reach a point and tolerance of reaching goal
goal_tol = 1e-3;
through_tol = 0.2;

% Size of the env't in (m) x = [-7.5m to 7.5m] and y = [-7.5m to 7.5m]
size_env = [15 15];

% Define percistent variables
persistent path size_map debugplot n;
% Where path contains through which the shortest path will passthrough
% size_map is the size of the map (pixel X pixel)
% debugplot to plot postition of the robot on the path while moving

% Initialize the robot control variables (returned by this function)
finish = false;
forwBackVel = 0;
leftRightVel = 0;
rotVel = 0;

% Initialization
if isempty(path)
    % Load image as a mtrix of pixel intensities
    map = double(imread(env_map)');

    % Size of the map
    size_map = size(map);

    % Start location on the map
    start_on_map = env2map(pos, size_map, size_env);

    % Goal location on the map
    goal_on_map = env2map(goal_pos, size_map, size_env);

    % Use wavefron algorithm to get a new map with wave connecting
    % start_on_map and goal_on_map
    map = wavefront(start_on_map, goal_on_map, map, debug);

    % Backpropagate the shortest path
    path = backpropagate(goal_on_map, map, debug);

    % Transform path from map to environment
    path = map2env(path, size_map, size_env);

    if debug
        hold on
        % Plot the moving robot as a green circle
        debugplot = plot(env2map(pos, size_map, size_env), 'go', MarkerSize=10);
        hold off
    end

elseif debug
    loc = env2map(pos, size_map, size_env);
    set(debugplot, 'Xdata', loc(1));
    set(debugplot, 'Ydata', loc(2));
end

% Follow the waypoints on the path
% Take the first cell on the path
diff_pos = path(end, :) - pos;
norm_diff = norm(diff_pos);

% Check if we have multiple waypoints on the path or this is the last cell
if size(path, 1) > 1
    % We have multiple points, we want to go through all of them with max
    % velocity magnetude, vel.
    vel = u_max;

    % check if a way point is reached
    if norm_diff < through_tol
        % If so remove the reached waypoint from the path
        path = path(1:end-1, :);

        % compute error for the next way point
        diff_pos = path(end, :) - pos;
        norm_diff = norm(diff_pos);
    end

else
    % This is the last waypoint on the path
    if norm_diff < goal_tol
        % Goal is reched
        % Stop the robot and clean exit
        leftRightVel = 0;
        forwBackVel = 0;
        rotVel = 0;
        finish = true;
        return;
    else
        % move with a velocity proportional to the error
        % Regulator output
        Pl = Kp_linear * norm_diff;
        % Limit Pl to get vel
        vel = max(min(Pl, u_max), -u_max);
    end

end

% Linear velocity computation
% Direction
direction = diff_pos/norm_diff;

% Velocity wrt global frame
u_global = vel * direction;

% Velocity wrt local frame
u_local = global2local(u_global, phi);

% Linear velocity commands
leftRightVel = u_local(1);
forwBackVel = u_local(2);

% Rotational velocity computation
% Make the robot face the direction of motion
diff_phi = angdiff(atan2(diff_pos(2), diff_pos(1)), phi - pi/2);

% Regulator output
Pr = Kp_angular * diff_phi;

% Limit Pr to get rotVel
rotVel = max(min(Pr, w_max), -w_max);

end


function on_map = env2map(on_env, size_map, size_env)
% Make a transformation from environment to map
on_map = round((on_env + size_env / 2) .* (size_map ./ size_env));
end

function on_env = map2env(on_map, size_map, size_env)
% Make a transformation from map to environment
on_env = on_map .* (size_env ./ size_map) - size_env / 2;
end

function u_local = global2local(u_global, phi)
% Transform velocity from global to local frame
% Rotation matrix from global to local frame
R = [cos(phi), -sin(phi);
    sin(phi), cos(phi)];

% Velocity in the local frame
u_local = R'*u_global';

end


function map = wavefront(start_on_map, goal_on_map, map, debug)
% This function computes a new map with wavefronts using wavefront
% algorith. It accepts start_on_map, goal_on_map, map, and a bool for
% debug plot

if debug
    % Initialize a figure
    fig = figure(1);
    colormap(fig, gray(2));
    ax = axes(fig);
    % clear axes
    cla(ax);
    % adjust limit to fit the map
    xlim(ax, [1 numrows(map)]);
    ylim(ax, [1 numcols(map)]);
    % plot
    imagesc(ax, 'CData', map');
    % label
    ylabel('y')
    xlabel('x')
    % pause to animate
    pause(0.1);
end

% map size
size_map = size(map);

% Define cost to travel between adjusent cells cost4 or cost8 can be
% defined. Here we are using cost8 as it gives better result.
cost8 = [sqrt(2) 1 sqrt(2);
    1       0       1;
    sqrt(2) 1 sqrt(1)];

% Enlarge obstacles to clear the robot from the wall 7x7 square is used
map = imerode(map, ones(7));

% Update debug plot
if debug
    imagesc(ax, 'CData', map');
    pause(0.1);
end

% Set walls as NaN and floor as Inf
map(map ~= 0) = Inf; % floor
map(map == 0) = NaN; % wall

% Put zero at the start
map(start_on_map(1), start_on_map(2)) = 0;

% Start the wave from the start_on_map
front = map == 0;  % start point of the wavefront
map_copy = map;    % Make a copy of the map

% Update debug plot
if debug
    numColors = 1.5 * mean(size_map);
    cmap = generate_color(numColors);
    colormap(fig, flip(cmap));
end

while true
    % Update the wavefront at each iteration
    if debug
        map_plot = map; map_plot(front) = 1;
        imagesc(ax, 'CData', map_plot', [0 1.5 * mean(size_map)]);
        pause(0.1);
    end

    % find index of front and clear front for reuse
    front_loc = find(front);
    front(:) = 0;

    % Loop over the front of the wave computing new cost values
    for k = 1:length(front_loc)
        % finding row and column of the current front location from idx
        [row, col] = ind2sub(size_map, front_loc(k));

        % 8 neighborhs of the front location
        neighbor_rows = max(row - 1, 1):min(row + 1, size_map(1));
        neighbor_cols = max(col - 1, 1):min(col + 1, size_map(2));

        % Compute cost values
        cost = map_copy(row, col) + cost8(neighbor_rows - row + 2, neighbor_cols - col + 2);

        % Apply the new map values if it reduces the cost
        map(neighbor_rows, neighbor_cols) = min(map(neighbor_rows, neighbor_cols), cost, 'includenan');

        % New fronts
        front(neighbor_rows, neighbor_cols) = max(map(neighbor_rows, neighbor_cols) == cost, front(neighbor_rows, neighbor_cols));
    end

    % Shouldn't repeat fronts
    front(front_loc) = 0;
    map_copy = map;

    % Check if goal is reached
    if ~isinf(map(goal_on_map(1), goal_on_map(2)))
        break;
    end
end
if debug
    imagesc(ax, 'CData', map', [0 1.5 * mean(size_map)]);
end

end

function path = backpropagate(goal_on_map, map, debug)
% Back propates the shortest path from goal to start using the cost of
% cells

% Extract the row and column indices of the goal cell from the input argument
row = goal_on_map(1);
col = goal_on_map(2);

% If debug is true, plot the goal cell as a black star on the map
if debug
    hold on
    plot(row, col, 'k*')
end

% Initialize the path as an empty matrix
path = [];

% While we haven't reached the start cell (which has a cost of 0)
while map(row, col) ~= 0
    % Add the current cell to the path
    path = [path; [row col]];
    % Mark the current cell as visited by setting its cost to infinity
    map(row, col) = Inf;

    % Find the indices of the 3x3 neighborhood around the current cell
    neighbor_rows = max(row - 1, 1):min(row + 1, numrows(map));
    neighbor_cols = max(col - 1, 1):min(col + 1, numcols(map));

    % Find the index of the cell in the neighborhood with the lowest cost
    [~, idx] = min(map(neighbor_rows, neighbor_cols), [], 'all', 'linear');

    % Convert the linear index to row and column indices
    [row_idx, col_idx] = ind2sub([3 3], idx);

    % Update the row and column indices to move to the lowest cost cell
    row = row + row_idx - 2;
    col = col + col_idx -2;

    % If we encounter a local minimum that is not the goal, raise an error
    if isinf(map(row, col))
        error('Could not reach the goal')
    end

end
% If debug is true, plot the shortest path as a red line on the map
if debug
    plot(path(:, 1), path(:, 2), 'r-', 'LineWidth', 2)
    hold off
end

end

function cmap = generate_color(numColors)
% This function returns a colormap for visualizing the wavefront algorithm in
% the debug plot. The colormap consists of a gradient from light blue to dark
% blue, with the last row representing walls.

% Define the RGB values for the light blue color
R0 = 173; G0 = 216; B0 = 230;

% Define the RGB values for the dark blue color
Rf = 23; Gf = 58; Bf = 85;

% Define the RGB values for the wall color (gray)
Rw = 128; Gw = 128; Bw = 128;

% Use linspace to create a set of evenly-spaced RGB values between the
% light blue color and the dark blue color. The number of colors in the colormap
% is specified by the input argument numColors.
cmap = [linspace(R0, Rf, numColors)' / 255, ...
    linspace(G0, Gf, numColors)' / 255, ...
    linspace(B0, Bf, numColors)' / 255];

% Set the last row of the colormap to gray to represent walls.
cmap(end, :) = [Rw, Gw, Bw] / 255;

end
