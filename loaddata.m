clear density dataname
dataname(1,:) = 'agpd1.out';
dataname(2,:) = 'agbx1.out';
dataname(3,:) = 'agbz1.out';
dataname(3,:) = 'agvx1.out';
dataname(4,:) = 'agde1.out';
numframes=50;
parfor i = 1:size(dataname,1);
density(i,:,:,:) = gasdatain(dataname(i,:),numframes,203);
end

