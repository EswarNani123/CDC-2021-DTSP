n_all = 3:13;
iter = 100;
rho = 10;
AngleAtFirstPoint = pi/6;
for i=1:length(n_all)
    n = n_all(i);
    OptimalTour = [];
    mincost = [];
    OptimalTourAngles = [];
    parfor j=1:iter 
        [OptimalTour(j,:),mincost(j),OptimalTourAngles(j,:)] = tsp_dp1(n,rho,AngleAtFirstPoint);
    end
    MeanDubinsTourCost(i) = mean(mincost);
end