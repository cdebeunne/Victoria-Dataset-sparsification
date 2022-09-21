%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%                      > x4 ---> x5
%                    /    ^
%                   /     |  
%        x1 ---> x2 --->  x3
import gtsam.*

graph = NonlinearFactorGraph;
priorNoise = noiseModel.Diagonal.Sigmas([0.3; 0.3; 0.1]);
graph.add(PriorFactorPose2(1, Pose2(0, 0, 0), priorNoise));

% Add odometry measurement 
noise = noiseModel.Diagonal.Sigmas([0.2; 0.2; 0.1]);
graph.add(BetweenFactorPose2(1, 2, Pose2(1, 0, pi/4), noise));
graph.add(BetweenFactorPose2(2, 3, Pose2(sqrt(0.5), -sqrt(0.5), -pi/4), noise));
graph.add(BetweenFactorPose2(2, 4, Pose2(sqrt(2), 0, -pi/4), noise));
graph.add(BetweenFactorPose2(3, 4, Pose2(0, 1, 0), noise));
graph.add(BetweenFactorPose2(4, 5, Pose2(1, 0, 0), noise));

% Set initial values
initial = Values;
initial.insert(1, Pose2(0.0, 0.0,  0.0));
initial.insert(2, Pose2(1.0, 0.0,  pi/4));
initial.insert(3, Pose2(2.0, 0.0, 0.0));
initial.insert(4, Pose2(2.0, 1.0, 0.0));
initial.insert(5, Pose2(3.0, 1.0, 0.0));

% Optimize using Levenberg-Marquardt optimization and get marginals
optimizer = LevenbergMarquardtOptimizer(graph, initial);
result = optimizer.optimizeSafely;
marginals = Marginals(graph, result);
plot2DTrajectory(result);



