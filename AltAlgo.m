function [DubinsTourCost, alpha, ElapsedTime] = AltAlgo(s,TourSeqDes,rho)
% rho = 1;
% n_all = 10:10:100;
% iter = 100;
ss = stateSpaceDubins;
ss.MinTurningRadius = rho;
DubinsTourCost = 0;
% dubConnObj = dubinsConnection;
% dubConnObj.MinTurningRadius = rho;
%for m=1:length(n_all)
%    n = n_all(m);
%    for i = 1:iter
%       DubinsTourCostAux = 0;
%       [TourSeqDes, ~, s] = TourSeqMinCost(n);
    tic;
        for j=1:(length(TourSeqDes)-1)
            if mod(j,2) == 1
                alpha(j) = atan2(s(TourSeqDes(j+1),2)-s(TourSeqDes(j),2),s(TourSeqDes(j+1),1)-s(TourSeqDes(j),1));
            else
                alpha(j) = alpha(j-1);
            end
        end
        alpha(end+1) = alpha(1);
        ElapsedTime = toc;
        for k=1:(length(TourSeqDes)-1)
            DubinsTourCost = DubinsTourCost+distance(ss, [s(TourSeqDes(k),:) alpha(k)], [s(TourSeqDes(k+1),:) alpha(k+1)]);
    %         [pathSegObj, ~] = connect(dubConnObj,[s(TourSeqDes(k),:) alpha(k)], [s(TourSeqDes(k+1),:) alpha(k+1)]);
    %         show(pathSegObj{1})
    %         hold on
        end
%        DubinsTourCost(i) = DubinsTourCostAux;
%    end
%   MeanDubinsTourCost(m) = mean(DubinsTourCost);
%end