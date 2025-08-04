

winsize = get(fig1,'Position');
%winsize(1:2) = [0 0];
A=moviein(25,fig1,winsize);
for i =1:25;
    imagep(:,:) = plas2(i,:,:);
    imageg(:,:) = gas2(i,:,:);
    subplot(h1)
    imagesc(imageg,clims1) 
    title('Normalised Gas Density','Fontsize',14,'Fontname','Arial')
    set(gca,'XTicklabel','')
    set(gca,'YTicklabel','')
    colorbar
    subplot(h2)
    imagesc(imagep,clims2) 
    title('Normalised Plasma Density','Fontsize',14,'Fontname','Arial')  
    set(gca,'XTicklabel','')
    set(gca,'YTicklabel','')
    colorbar
    A(:,i)=getframe(fig1,winsize);
end

