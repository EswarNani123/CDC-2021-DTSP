function [DubinsTourCost, alpha, ElapsedTime] = SVA(s,TourSeqDes,rho,AngleAtFirstPoint)
% rho = 0.1;
% n_all = 3:13;
% iter = 100;
ss = stateSpaceDubins;
ss.MinTurningRadius = rho;
dubConnObj = dubinsConnection;
dubConnObj.MinTurningRadius = rho;
%AngleAtFirstPoint = pi/6;
alpha(1) = AngleAtFirstPoint;
%for m=1:length(n_all)
n = length(s);
DubinsTourCost = 0;
%   for i = 1:iter
%       DubinsTourCostAux = 0;
%        [TourSeqDes, ~, s] = TourSeqMinCost(n);
        tic;
        for j=2:(length(TourSeqDes)-1)
            alpha(j) = final_heading([s(TourSeqDes(j-1),:) alpha(j-1)],s(TourSeqDes(j),:),rho);
        end
        alpha(end+1) = final_heading([s(TourSeqDes(end-1),:) alpha(n)],s(TourSeqDes(end),:),rho);
        ElapsedTime = toc;
        for k=1:(length(TourSeqDes)-1)
            DubinsTourCost = DubinsTourCost+distance(ss, [s(TourSeqDes(k),:) alpha(k)], [s(TourSeqDes(k+1),:) alpha(k+1)]);
%             [pathSegObj, ~] = connect(dubConnObj,[s(TourSeqDes(k),:) alpha(k)], [s(TourSeqDes(k+1),:) alpha(k+1)]);
%             show(pathSegObj{1})
%             hold on
        end
%        DubinsTourCost(i) = DubinsTourCostAux;
    %end
%   MeanDubinsTourCost(m) = mean(DubinsTourCost);
%end