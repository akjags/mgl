% mglDoubleBars.m
%
%        usage: mglDoubleBars(numBlocks)
%           by: akshay jagadeesh
%         date: 12/15/2016
%      purpose: Displays a retinotopy stimulus with two flashing checkerboard bars.
%
%

function mglDoubleBars(numBlocks)

global myScreen xLimit yLimit;
myScreen.background = 'gray';
myScreen.transparentBackground = false;
myScreen.setVolume = false;
myScreen = initScreen(myScreen);
xLimit = myScreen.imageWidth / 2;
yLimit = myScreen.imageHeight / 2;

% Draw fixation cross
mglFixationCross();

%frameNum = 4;
bar1Contrast = 0.3;
bar2Contrast = 0.6;

% Bar Directions
deg = 0:45:315;
dirs = [1 4; 1 8; 6 4; 6 8; 2 3; 2 7; 6 3; 6 7; 3 4; 3 8; 7 4; 7 8];
%dirs = dirs(randperm(length(dirs)),:);
for k = 1:8
dir1 = dirs(k,1);
dir2 = dirs(k,2);
disp(sprintf('Bar 1: (%i) %i degrees; Bar 2: (%i) %i degrees', dir1, deg(dir1), dir2, deg(dir2)));
  for frameNum = 1:20
  % Draw stencil for bar #1
  mglStencilCreateBegin(1);
  c1 = drawBars(dir1, frameNum);
  mglStencilCreateEnd;
  mglClearScreen;

  % Draw stencil for bar #2
  mglStencilCreateBegin(2);
  c2 = drawBars(dir2, frameNum);
  mglStencilCreateEnd;
  mglClearScreen;

  % Calculate overlap of bars
  %
  mglStencilCreateBegin(3);
  if dir1 == 1 || dir1 == 5
    calcDrawOverlap(c1,c2,1);
  elseif dir2 == 1 || dir2 == 5
    calcDrawOverlap(c1,c2,2);
  else
    calcDrawOverlap(c1, c2);
  end
  mglStencilCreateEnd;
  mglClearScreen;

    for i = 1:5
    % Draw checkerboard, stenciled by #1
    mglStencilSelect(1);
    drawCheckerboard(bar1Contrast, i);

    mglStencilSelect(2);
    drawCheckerboard(bar2Contrast, i);

    mglStencilSelect(3);
    drawCheckerboard(bar1Contrast + bar2Contrast, i);

    mglStencilSelect(0);
    mglFlush;
    mglClearScreen;
    end
  end
end


return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             Helper Methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ~~~ calcDrawOverlap ~~~ %%%%%%
% Calculates the coordinates of the area of overlap between two bars, and draws it using mglPolygon
% Arguments: 
%        bar1, bar2 --> coordinates arrays in the form [x1 y1; x2 y2; x3 y3; x4 y4]
function calcDrawOverlap(bar1, bar2, vert)

if ieNotDefined('vert')
  vert = 0;
end

if vert ~= 1
  % for each bar, first calculate the line (y=mx+b) between points 1 and 2, and points 3 and 4 (long sides)
  [rB1(1), rB1(2)] = calculateLine(bar1(1,:), bar1(2,:));
  [fB1(1), fB1(2)] = calculateLine(bar1(3,:), bar1(4,:));
end

if vert ~= 2
  [rB2(1), rB2(2)] = calculateLine(bar2(1,:), bar2(2,:));
  [fB2(1), fB2(2)] = calculateLine(bar2(3,:), bar2(4,:));
end

if vert == 1
  x1 = bar1(1,1); x2 = bar1(3,1);
  p1(1) = x1; p1(2) = rB2(1)*x1 + rB2(2);
  p2(1) = x2; p2(2) = rB2(1)*x2 + rB2(2);
  p3(1) = x2; p3(2) = fB2(1)*x2 + fB2(2);
  p4(1) = x1; p4(2) = fB2(1)*x1 + fB2(2);
elseif vert == 2
  x1 = bar2(1,1); x2 = bar2(1,1);
  p1(1) = x1; p1(2) = rB1(1)*x1 + rB1(2);
  p2(1) = x2; p2(2) = rB1(1)*x2 + rB1(2);
  p3(1) = x2; p3(2) = fB1(1)*x2 + fB1(2);
  p4(1) = x1; p4(2) = fB1(1)*x1 + fB1(2);
else
  [p1(1), p1(2)] = calculateIntersect(rB1, rB2);
  [p2(1), p2(2)] = calculateIntersect(rB1, fB2);
  [p3(1), p3(2)] = calculateIntersect(fB1, fB2);
  [p4(1), p4(2)] = calculateIntersect(fB1, rB2);
end

intersect = [p1; p2; p3; p4];
mglPolygon(intersect(:,1)', intersect(:,2)', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ~~~ calculateIntersect ~~~ %%%%%%
% Calculates the point of intersection given two lines
%   Arguments: 
%        l1, l2 --> lines in the form [slope, intersect]
function [intersectX,intersectY] = calculateIntersect(l1, l2)
intersectX = (l2(2) - l1(2)) / (l1(1) - l2(1));
intersectY = (l1(1)*l2(2) - l2(1)*l1(2)) / (l1(1) - l2(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ~~~ calculateLine ~~~ %%%%%%%%
% Calculates a line between two points
%   Arguments: 
%        p1, p2 --> points in the form [x,y]
function [slope,intercept] = calculateLine(p1, p2)
slope = (p2(2) - p1(2)) / (p2(1) - p1(1));
intercept = p2(2) - slope*p2(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ~~~ drawCheckerBoard ~~~ %%%%%%
% Draws a checkerboard to fill the screen using mglQuad
%   Starting in the top left, draws alternating black and white squares.
%   Then imposes a gray (0.50) semi-transparent window
%   Arguments: 
%        contrast --> 1-alpha : specifies transparency of window, i.e. contrast of checkerboard. 
%        index --> 0 or 1 specifies which color to start with (for flashing effect)
%
function drawCheckerboard(contrast, index)

global myScreen xLimit yLimit;
numWide = 41; numTall = 29;

sqWid = xLimit/20;
color = mod(index,2);
for i = 0:numWide
  for j = 0:numTall
    xCoords = [-xLimit + i*sqWid; -xLimit + (i+1)*sqWid; -xLimit + (i+1)*sqWid; -xLimit + i*sqWid];
    yCoords = [-yLimit + (j+1)*sqWid; -yLimit + (j+1)*sqWid; -yLimit + j*sqWid; -yLimit + j*sqWid];
    mglPolygon(xCoords', yCoords', color);
    color = mod(color+1,2);
  end
  color = mod(color+1,2); 
end
% superimpose a semitransparent polygon across the whole screen with alpha = 1-contrast
mglPolygon([-xLimit, xLimit, xLimit, -xLimit], [yLimit, yLimit, -yLimit, -yLimit], [0.5 0.5 0.5 1-contrast]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ~~~ drawBars ~~~ %%%%%%%%%%%%%
% Draw a single bar with width 2, according to the direction and frameNum
%   Arguments:
%        barDirection --> Numbers 1 through 8 map to directions from 0 degrees (vertical bar 
%                         moving to right) to 315 (bar with slope -1 moving from bottom left to top right)
%        frameNum --> Specifies how far the bar has swept across the screen from 0 (initial position) to 20 (final position)
%
function coords = drawBars(barDirection, frameNum)
if ieNotDefined('frameNum')
  frameNum = 0;
end
% Initialize bars 10 steps behind center
init = 10;

global myScreen xLimit yLimit;
switch barDirection
  case 1 % 0 degrees: vertical bar moves left to right
    xCorners = [-xLimit+frameNum, -xLimit+frameNum, -xLimit+2+frameNum, -xLimit+2+frameNum]; 
    yCorners = [yLimit, -yLimit, -yLimit, yLimit];
  case 2 % 45: bar moves top left to bottom right
    xCorners = [xLimit - sqrt(2)-init+frameNum, -xLimit-init+frameNum, -xLimit+sqrt(2)-init+frameNum, xLimit-init+frameNum];
    yCorners = [yLimit+init-frameNum, -yLimit+sqrt(2)+init-frameNum, -yLimit+init-frameNum, yLimit-sqrt(2)+init-frameNum];
  case 3
    xCorners = [-xLimit, xLimit, xLimit, -xLimit];
    yCorners = [yLimit-frameNum, yLimit-frameNum, yLimit-2-frameNum, yLimit-2-frameNum];
  case 4 % 135: bar moves top right to bottom left
    xCorners = [-xLimit+sqrt(2)+init-frameNum, xLimit+init-frameNum, xLimit-sqrt(2)+init-frameNum, -xLimit+init-frameNum];
    yCorners = [yLimit+init-frameNum, -yLimit+sqrt(2)+init-frameNum, -yLimit+init-frameNum, yLimit-sqrt(2)+init-frameNum];
  case 5
    xCorners = [xLimit-frameNum, xLimit-frameNum, xLimit-2-frameNum, xLimit-2-frameNum];
    yCorners = [yLimit, -yLimit, -yLimit, yLimit];
  case 6
    xCorners = [xLimit+init-frameNum, -xLimit+sqrt(2)+init-frameNum, -xLimit+init-frameNum, xLimit-sqrt(2)+init-frameNum];
    yCorners = [yLimit-sqrt(2)-init+frameNum, -yLimit-init+frameNum, -yLimit+sqrt(2)-init+frameNum, yLimit-init+frameNum];
  case 7
    xCorners = [-xLimit, xLimit, xLimit, -xLimit]; 
    yCorners = [-yLimit+frameNum, -yLimit+frameNum, -yLimit+2+frameNum, -yLimit+2+frameNum];
  case 8
    xCorners = [-xLimit-init+frameNum, xLimit-sqrt(2)-init+frameNum, xLimit-init+frameNum, -xLimit+sqrt(2)-init+frameNum];
    yCorners = [yLimit-sqrt(2)-init+frameNum, -yLimit-init+frameNum, -yLimit+sqrt(2)-init+frameNum, yLimit-init+frameNum];
end
mglPolygon(xCorners, yCorners, 0);
coords = [xCorners; yCorners]';
