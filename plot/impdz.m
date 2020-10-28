function tt1g=imtz(t,clims)
    L=2*pi;
    anom=true;
    cmap=redblue();
    %clims=[-0.5 0.5];
    fid = fopen("coord");
    coord = textscan(fid, "%s %d %s %f");
    fclose(fid);
    nz=length(coord{:,2});
    nz=floor((nz-1)/2)+1;
    coord=coord{1:nz,4};
    nz=nz-2;
    filename=sprintf('p%07d.zzz',t);
    tt=load(filename); 
	n=length(tt); n=n/nz;
	zz=coord(2:(nz+1));
    xx=(0:(n-1))/n*L;
	tt=reshape(tt,n,nz);
    if (anom)
	    tt1=tt-ones(n,1)*mean(tt);
    else
        tt1=tt;
    end
    tt1g=(tt1(:,2:nz)-tt1(:,1:(nz-1)))./(ones(n,1)*(zz(2:nz)-zz(1:(nz-1)))');
    pcolor(xx,zz(2:(nz-1))+0.5 , tt1g(:,1:(nz-2))'); shading flat; hc=colorbar;
    ylabel(hc,'T')
    caxis(clims);
    colormap(cmap);
    set(gca,'YDir','normal','FontSize',12);
    set(gca,'XTick',[0 pi/2 pi 1.5*pi (pi*2)-(2*pi/n)])
    set(gca,'XTickLabel',{'0', '\pi/2', '\pi', '3/2\pi', '2\pi'})
    xlabel('x');
    ylabel('z');
    title(sprintf("step=%d",t));
    set(gca,'Fontsize',16)
end
    
