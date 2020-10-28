function imtz(t,clims)
    L=2*pi;
    anom=false;
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
    pp=load(filename); 
    filename=sprintf('t%07d.zzz',t);
    tt=load(filename); 
	n=length(pp); n=n/nz;
	zz=coord(2:(nz+1));
    xx=(0:(n-1))/n*L;
	pp=reshape(pp,n,nz);
	tt=reshape(tt,n,nz);

    if (anom)
	    pp1=pp-ones(n,1)*mean(pp);
	    tt1=tt-ones(n,1)*mean(tt);
    else
        pp1=pp;
        tt1=tt;
    end
    pp1g=(pp1(:,2:(nz-1))-pp1(:,1:(nz-2)))./(ones(n,1)*(zz(2:(nz-1))-zz(1:(nz-2)))');
    tt1g=1e8*tt(:,2:(nz-1));
    pcolor(xx(2:(n-1)),zz(2:(nz-1))+0.5 , pp1g(2:(n-1),:)'-tt1g(2:(n-1),:)'); shading flat; hc=colorbar;
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
    
