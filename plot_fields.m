function [fig, ax] = plot_fields(X, Y, u, v, ux, uy, vx, vy, EdMag, EdPhi, fntSz, linSz)
%PLOT_FIELDS Plot order, velocity, strain-rate fields, and stress fields.
%   X and Y are matching grid matrices; fields follow X(:) ordering.
%   Returns the figure and axes for adjusting color limits before export.
if nargin < 11, fntSz = 20; end
if nargin < 12, linSz = 1; end

%% Plot settings
coarse = 7;                    % Grid stride for vector overlays
boundaryWidth = 3;
theta = linspace(0, 2*pi, 6290);
circleX = cos(theta);
circleY = sin(theta);

%% Prepare fields on the plotting grid
gridSize = size(X);
u = reshape(u, gridSize);
v = reshape(v, gridSize);
ux = reshape(ux, gridSize);
uy = reshape(uy, gridSize);
vx = reshape(vx, gridSize);
vy = reshape(vy, gridSize);
strainMagnitude = reshape(EdMag, gridSize);
strainAngle = reshape(EdPhi, gridSize);

[orderMagnitude, directorAngle] = orient_order(X, Y);
speed = hypot(u, v);
isotropicStrainRate = 0.5 * (ux + vy);

beta1 = isotropicStrainRate;
beta2 = 0.5 * (ux - vy);
beta3 = 0.5 * (uy + vx);
[isotropicStress, stressDevXX, stressDevXY] = compute_stress(beta1, beta2, beta3, X, Y);
[material, kappa, prestress] = get_parameters;
meanStress = isotropicStress + prestress(1);

% Deviatoric stress is [stressDevXX, stressDevXY; stressDevXY, -stressDevXX].
% Its positive eigenvalue gives the magnitude; the eigenvector gives the axis.
devStressMagnitude = hypot(stressDevXX, stressDevXY);
devStressAngle = 0.5 * atan2(stressDevXY, stressDevXX);
devStressAngle(devStressMagnitude == 0) = NaN; % Axis is undefined at zero stress

%% Define panels in display order
scalarFields = {orderMagnitude, speed, isotropicStrainRate, strainMagnitude, meanStress, devStressMagnitude};
vectorX = {cos(directorAngle), u, u, strainMagnitude .* cos(strainAngle), u, ...
    devStressMagnitude .* cos(devStressAngle)};
vectorY = {sin(directorAngle), v, v, strainMagnitude .* sin(strainAngle), v, ...
    devStressMagnitude .* sin(devStressAngle)};
panelTitles = {'$Q(\mathbf{x})$', '$\mathbf{u}(\mathbf{x})$', ...
    '$\mathrm{Tr}[\nabla\mathbf{u}]/2$', '$e_{\mathrm{d}},\ \mathbf{p}$', ...
    '$\mathrm{Tr}[\sigma]/2$', '$\sigma_{\mathrm{d}},\ \mathbf{p}_{\sigma}$'};
vectorColors = {'m', 'k', 'k', 'm', 'k', 'm'};
arrowHeads = {'off', 'on', 'on', 'off', 'on', 'off'};

% Use the same unit-disk mask for contours and vector overlays.
outside = hypot(X, Y) > 1;
rows = 1:coarse:size(X, 1);
cols = 1:coarse:size(X, 2);
coarseX = X(rows, cols);
coarseY = Y(rows, cols);

%% Draw panels with shared formatting
fig = figure('Units', 'normalized', 'Color', 'w', ...
    'Position', [0.0250, 0.3236, 0.9497, 0.3866]);
panelCount = numel(scalarFields);
layout = tiledlayout(fig, 1, panelCount, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(1, panelCount);

for panel = 1:panelCount
    field = scalarFields{panel};
    dx = vectorX{panel};
    dy = vectorY{panel};
    field(outside) = NaN;
    dx(outside) = NaN;
    dy(outside) = NaN;

    ax(panel) = nexttile(layout);
    contourf(ax(panel), X, Y, field, 'EdgeColor', 'none');
    hold(ax(panel), 'on');
    quiver(ax(panel), coarseX, coarseY, dx(rows, cols), dy(rows, cols), ...
        'Color', vectorColors{panel}, 'LineWidth', linSz, ...
        'ShowArrowHead', arrowHeads{panel});
    plot(ax(panel), circleX, circleY, 'k', 'LineWidth', boundaryWidth);
    hold(ax(panel), 'off');
    format_panel(ax(panel), panelTitles{panel}, fntSz);
end

% Double backslashes because sprintf interprets escape sequences.
parameterTitle = sprintf([ ...
    '$\\nu = %.1f$, $\\lambda^o = %.1f$, ' ...
    '$\\lambda^s = %.1f$, $\\delta\\nu = %.1f$, ' ...
    '$m_0 = %.1f$'], ...
    material(2), material(3), material(4), material(5), ...
    prestress(1));
title(layout, parameterTitle, 'Interpreter', 'latex', ...
    'FontSize', 24, 'Color', 'k');
end

function format_panel(ax, panelTitle, fontSize)
% Keep limits, typography, and colorbar placement consistent.
axis(ax, 'equal');
axis(ax, 'off');
xlim(ax, [-1.05, 1.05]);
ylim(ax, [-1.05, 1.05]);
set(ax, 'FontSize', fontSize, 'Color', 'w');
colormap(ax, parula(256));
title(ax, panelTitle, 'Interpreter', 'latex', 'FontSize', fontSize, 'Color', 'k');
cb = colorbar(ax, 'southoutside');
cb.FontSize = fontSize;
cb.LineWidth = 3;
cb.Color = 'k';
end
