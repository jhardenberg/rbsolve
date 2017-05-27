	function im(t,nz,sc)
	%clf
    %t
	load coord1
	filename=sprintf('t%07d.zzz',t);
        tt=load(filename); 
	n=length(tt); n=n/nz;
	zz=coord1(2:(nz+1),2);
        xx=(0:(n-1))/n*2*pi;
	tt=reshape(tt,n,nz);
	tt1=tt;%-ones(n,1)*mean(tt);
	%  subplot(2,1,nn)
	if(length(sc)==0)
	pcolor(xx,zz , tt1'); shading flat; colorbar;
    else
	pcolor(xx,zz, tt1'); shading flat;colorbar;
        end
	set(gca,'YDir','normal','FontSize',12);
	xlabel('x');
	ylabel('z');
