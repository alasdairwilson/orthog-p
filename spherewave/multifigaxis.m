for i = 1 : 4;
subaxis(2,2,i, 'Spacing', 0.03, 'Padding', 0.01, 'Margin', 0.1);
imagesc(A2(:,:,i)*10,[-0,0.0001])
axis off
axis equal
axis on
xlim([1,200])
xlabel('x','FontSize',18,'Rotation',0,'FontName','Times New Roman')
ylabel('y','FontSize',18,'Rotation',0,'FontName','Times New Roman')
title(timetitle(i,:),'Fontsize',20,'FontName','Times New Roman')
colorbar
end

p=mtit(thefig,'Evolution of Fresh Plasma Created by AI','Fontsize',22,'FontName','Arial');