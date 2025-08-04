fig1 = figure(1)
subplot(3,3,1)
image(:,:) = gas(6,:,:);
imagesc(image,[0.997,1.003])
axis square
originalsize1 = get(gca, 'Position');

subplot(3,3,2)
image(:,:) = gas(10,:,:);
imagesc(image,[0.997,1.003])
axis square
subplot(3,3,3)
image(:,:) = gas(14,:,:);
imagesc(image,[0.997,1.003])
axis square
colorbar('YTickLabel',{'0.97','1.00','1.03'})
thisposition = get(gca, 'Position');
thisposition(3:4) = originalsize1(3:4);
set(gca,'Position',thisposition)


subplot(3,3,4)
image(:,:) = gasi(6,:,:);
imagesc(image,[0.997,1.003])
axis square
subplot(3,3,5)
image(:,:) = gasi(10,:,:);
imagesc(image,[0.997,1.003])
axis square
subplot(3,3,6)
image(:,:) = gasi(14,:,:);
imagesc(image,[0.997,1.003])
axis square
colorbar('YTickLabel',{'0.97','1.00','1.03'})
thisposition = get(gca, 'Position');
thisposition(3:4) = originalsize1(3:4);
set(gca,'Position',thisposition)

subplot(3,3,7)
image(:,:) = diff(6,:,:);
imagesc(image,[-0.004,0.004])
axis square
subplot(3,3,8)
image(:,:) = diff(10,:,:);
imagesc(image,[-0.004,0.004])
axis square
subplot(3,3,9)
image(:,:) = diff(14,:,:);
imagesc(image,[-0.004,0.004])
axis square
colorbar('YTickLabel',{'-0.004','-0.002','0.000','0.002','0.004'})
thisposition = get(gca, 'Position');
thisposition(3:4) = originalsize1(3:4);
set(gca,'Position',thisposition)