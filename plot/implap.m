function plap=impz(t, clims)
    L=2*pi;
    anom=false;
    cmap=redblue();
    fid = fopen("coord");
    coord = textscan(fid, "%s %d %s %f");
    fclose(fid);
    nz=length(coord{:,2});
    nz=floor((nz-1)/2)+1;
    zz=coord{1,4}((nz+1):(end));
    nz=nz-3;
    filename=sprintf('p%07d.zzz',t);
    tt=load(filename); 
    n=length(tt); n=n/(nz+1);
    xx=(0:(n-1))/n*L;
    tt=reshape(tt,n,nz+1); tt=tt(:,1:nz);
    if (anom)
	    tt1=tt-ones(n,1)*mean(tt);
    else
        tt1=tt;
    end
    dx=xx(2)-xx(1);
    tt1x=(tt1(3:n,2:(nz-1))+tt1(1:(n-2),2:(nz-1))-2*tt1(2:(n-1),2:(nz-1)))./(dx*dx);
    
    %dza=(ones(n-2,1)*(zz(3:nz)-zz(2:(nz-1)))');
    %dzb=(ones(n-2,1)*(zz(2:(nz-1))-zz(1:(nz-2)))');
    %tt1a=2*(tt1(2:(n-1),3:nz)-tt1(2:(n-1),2:(nz-1)))./(dza.*(dza+dzb));
    %tt1b=2*(tt1(2:(n-1),1:(nz-2))-tt1(2:(n-1),2:(nz-1)))./(dzb.*(dza+dzb))

    yp=zz;

    j=2:(nz+1);
    qp=2.d0./((yp(j+1)-yp(j-1)).*(yp(j)-yp(j-1))); qp=ones((n-2),1)*qp';
    rp=2.d0./((yp(j+1)-yp(j)).*(yp(j)-yp(j-1))); rp=ones((n-2),1)*rp';
    sp=2.d0./((yp(j+1)-yp(j-1)).*(yp(j+1)-yp(j))); sp=ones((n-2),1)*sp';
    tt1a=tt1(2:(n-1),1:(nz-2)).*qp(:,2:(end-1))- tt1(2:(n-1),2:(nz-1)).*rp(:,2:(end-1))+ tt1(2:(n-1),3:(nz)).*sp(:,2:(end-1));
    plap=tt1a+tt1x;
   % plap=tt1x;

    pcolor(xx(2:(n-1)),zz(2:(nz-1))+0.5 , plap'); shading flat; hc=colorbar;
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
    
