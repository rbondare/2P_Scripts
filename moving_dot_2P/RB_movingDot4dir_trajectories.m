function plan = RB_movingDot4dir_trajectories(sp, canvas_wh, full_screen_az, ifi)
% RB_MOVINGDOT4DIR_TRAJECTORIES  Trajectory table of the 'moving_dot_4dir' stimulus.
%
%   plan = RB_movingDot4dir_trajectories(specParams, [width height], full_screen_az, ifi)
%
% Geometry and timing of the dot sweeps used by present_stimulus_addnew
% (case 'moving_dot_4dir'), ported from moving_dot_4dir_opto_oncheckers.m.
% A dot starts on one wall of a rectangular area, rests, crosses the area in a
% straight line to the opposite wall, and rests again. Starting positions lie one
% dot radius apart along the walls, and every position is used in the two
% directions along its axis, so the area is covered in 4 cardinal directions.
% One trajectory = one subtrial.
%
% specParams fields used
%   dot_speed_deg_s, dot_radius_deg   dot speed and radius
%   dot_start_t, dot_stop_t           s of rest before / after the motion
%                                     (dot_start_t is forced to >= 0.1 s so that
%                                     directions can be separated afterwards)
%   visual_area_10_20                 [bottom top left right] on a 10 x 20 grid
%                                     of the canvas, e.g. [0 5 10 20]
%
% plan fields (pixels are canvas pixels, x = column, y = row)
%   traj            [n x 6]  x0 y0 vx vy x1 y1, v in px per frame
%   nsteps          [n x 1]  frames spent moving
%   nframes         [n x 1]  frames of the whole trajectory (rest + motion + rest)
%   dir_deg         [n x 1]  motion direction on the canvas: 0 right, 90 down,
%                            180 left, 270 up (y grows downwards)
%   tstart_frames, tstop_frames, radius_px, step_px, area_px ([ymin ymax xmin xmax]), n
%
% Position in frame k (0-based) of trajectory i, as in the original script:
%   traj(i,1:2) + min(max(k-tstart_frames,0), nsteps(i)) * traj(i,3:4)

w = canvas_wh(1);
h = canvas_wh(2);
px_per_deg = w/full_screen_az;
rpx = max(1, round(sp.dot_radius_deg*px_per_deg)); %dot radius in px
dx = sp.dot_speed_deg_s*px_per_deg*ifi;            %dot speed in px per frame
delta_pos = rpx;                                   %neighbouring trajectories are a radius apart

%area to cover in px, stretched to a whole number of trajectory spacings
area = [sp.visual_area_10_20(1:2)/10*h, sp.visual_area_10_20(3:4)/20*w];
area(2) = area(1)+ceil(diff(area(1:2))/delta_pos)*delta_pos;
area(4) = area(3)+ceil(diff(area(3:4))/delta_pos)*delta_pos;
ys = area(1):delta_pos:area(2);
xs = area(3):delta_pos:area(4);
nhz = numel(ys); %horizontal trajectories
nvr = numel(xs); %vertical trajectories

%x0 y0 vx vy x1 y1
traj = [repmat(area(3),nhz,1) ys(:) repmat([ dx 0],nhz,1) repmat(area(4),nhz,1) ys(:);   %left wall -> right
        repmat(area(4),nhz,1) ys(:) repmat([-dx 0],nhz,1) repmat(area(3),nhz,1) ys(:);   %right wall -> left
        xs(:) repmat(area(2),nvr,1) repmat([0 -dx],nvr,1) xs(:) repmat(area(1),nvr,1);   %y=area(2) wall -> y=area(1)
        xs(:) repmat(area(1),nvr,1) repmat([0  dx],nvr,1) xs(:) repmat(area(2),nvr,1)];  %y=area(1) wall -> y=area(2)

tstart_frames = round(max(0.1, sp.dot_start_t)/ifi);
tstop_frames = round(sp.dot_stop_t/ifi);
nsteps = ceil((abs(traj(:,5)-traj(:,1))+abs(traj(:,6)-traj(:,2)))/dx);

plan.traj = traj;
plan.nsteps = nsteps;
plan.nframes = tstart_frames+tstop_frames+nsteps;
plan.dir_deg = mod(atan2d(traj(:,4),traj(:,3)),360);
plan.tstart_frames = tstart_frames;
plan.tstop_frames = tstop_frames;
plan.radius_px = rpx;
plan.step_px = dx;
plan.area_px = area;
plan.n = size(traj,1);
end
