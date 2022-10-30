function [ky,kz] = ROCK(demosig,nav_interval,k,matrixsize,lambda)
% matrixsize = [Nkz,Nky]

Nkz = matrixsize(1); Nky =  matrixsize(2);
golden_ratio = (sqrt(5)+1)/2;
golden_angle = pi/golden_ratio;


numRings = nav_interval-1;

kcenter = floor(matrixsize/2)+1;

% kz^^^kz  ky>>>ky
kykz = zeros(matrixsize);


ky = zeros(size(demosig));
kz = zeros(size(demosig));

[K1,K2] = meshgrid(1:5);
K = [K1(:) K2(:)];
sigma = eye(2);
kernel = mvnpdf(K,[3,3],sigma);
kernel = reshape(kernel,5,5);

[Y,Z] = meshgrid(1:Nky,1:Nkz);
for i = 1:numRings
    
    %ring elements
    distToCenter = sqrt(((Z-kcenter(1))/(0.5*Nkz)).^2+((Y-kcenter(2))/(0.5*Nky)).^2);
    ringele = distToCenter>(i-1)/numRings & distToCenter<=(i)/numRings;
    
    ringind = find(ringele);
    
    [zind,yind] = ind2sub(matrixsize,ringind);
    
    rings{i} = [zind,yind];
    
end

densMtx = zeros(matrixsize);
for i = demosig
    if mod(i-1,nav_interval) == 0
        ky(i) = kcenter(2);
        kz(i) = kcenter(1);
    else
        
        armNum = floor((i-1)/nav_interval);
        PHI = mod(armNum*golden_angle,2*pi);
        armOrder = mod(i-2,nav_interval)+1; % 1 to nav_interval-1
        
        %s spiral in, reverse orderq
        armOrder = (nav_interval)-armOrder; % nav_interval-1 to 1
%         r = armOrder./nav_interval;
        
        ringOfi = rings{armOrder};
        
        zind = ringOfi(:,1);
        yind = ringOfi(:,2);
        
        C = Inf;
        
        
        

        for j = 1:length(yind)
            complexVal = complex((yind(j)-kcenter(2))/(0.5*Nky),(zind(j)-kcenter(1))/(0.5*Nkz));
            thetayz = angle(complexVal);
            r = abs(complexVal);
            phiyz = mod(thetayz-r*k,2*pi);
            
            tempC = min(abs(phiyz-PHI),2*pi-abs(phiyz-PHI))+lambda*densMtx(zind(j),yind(j));
            
            if tempC<C
               C = tempC;
               ymin = yind(j);
               zmin = zind(j);
            end
        end
        
        
        ky(i) = ymin;
        kz(i) = zmin;
        
        oneMtx = zeros(matrixsize);
        oneMtx(zmin,ymin) = 1;
        
        densMtx = densMtx+conv2(oneMtx,kernel,'same');
    end
    
    
end
end

function density = D(densMtx,yi,zi,radius)
    [X,Y] = meshgrid(-radius:radius);
    
    circle = sqrt(X.^2+Y.^2)<=radius;
    
    surroundInd = find(circle);
%     [xSurInd,ySurInd] = ind2sub(size(X),surroundInd);
    density = 0;
    for i = 1:length(surroundInd)
        density = density+densMtx(zi+Y(surroundInd(i)),yi+X(surroundInd(i)));
    end

end