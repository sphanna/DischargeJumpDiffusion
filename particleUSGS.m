T = readtable('USGS2015-2020.csv');
discharge = T.x293221_00060_00003;
time = T.datetime;

pf = particleFilter(@particleUSGS_Propagator,@particleLikelihoodFcn);
initialize(pf, 5000, [200;2;0.7;0.1;0.04;0.025;0.06;200], 0.1*eye(8));

l = size(discharge,1);
N = 500;
result = zeros([N,8]);
dischargeSet = discharge(1:N);
for k = 1:N
    result(k,:) = correct(pf,discharge(k));
    predict(pf);
end

result(N,:)

begin = 1; en = N;
xlabel('Time (Days)') 
ylabel('Discharge (ft^3 / s)') 
title('Discharge over Time: Filtered Results');
plot(result(begin:en,1),'red'); hold on
plot(discharge(begin:en),'blue'); hold off

err = (result(:,1) - discharge)
xlabel('Time (Days)') 
ylabel('Discharge (ft^3 / s)') 
title('Residuals');
plot(err)

xlabel('Time (Days)') 
ylabel('Saturation Percentage') 
title('Saturation Percentage Parameter (p)');
plot(result(:,2))

xlabel('Time (Days)') 
ylabel('Rate') 
title('Precipitation Occurrance Rate (lambda)');
plot(result(:,5))

xlabel('Time (Days)') 
ylabel('Precipitation Magnitude') 
title('Precipitation Magnitude Parameter (M)');
plot(result(:,8))

function likelihood = particleLikelihoodFcn(particles, measurement)

    yHat = particles(1,:);
    
    e = bsxfun(@minus, yHat, measurement(:)'); % Error
    numberOfMeasurements = 1;
    mu = 0; % Mean
    Sigma = eye(numberOfMeasurements); % Variance
    measurementErrorProd = dot((e-mu), Sigma \ (e-mu), 1);
    c = 1/sqrt((2*pi)^numberOfMeasurements * det(Sigma));
    likelihood = c * exp(-0.5 * measurementErrorProd);
end

function xNext = particleUSGS_Propagator(x)
    x = abs(x);
    xNext = x;
    [numParams,numParticles] = size(xNext);
    jump = zeros([1,numParticles]);
    for j = 1:numParticles
        jump(j) = sum(exprnd(x(8,j),[1,poissrnd(x(7,j))]));
    end
    
    xNext(1,:) = x(2,:) - bsxfun(@times,x(3,:),x(1,:)) + bsxfun(@times,x(4,:),randn([1,numParticles])) + jump;
    xNext(2,:) = bsxfun(@times,x(5,:),x(1,:)) - bsxfun(@times,x(6,:),x(2,:));
    
    %parameter error
    xNext = xNext + bsxfun(@times,[0.1; 0.1; 1e-1; 1e-1; 1e-2; 1e-2; 1e-2; 1e2],randn(size(xNext)));
    xNext = abs(xNext);
end
