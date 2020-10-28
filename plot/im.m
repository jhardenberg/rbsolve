	function im(t,sc)

	filename=sprintf('t%07d.dat',t);
        tt=load(filename); 
	n=length(tt); n=sqrt(n);
	tt=reshape(tt,n,n);
%	tt=tt-mean(mean(tt));
	imagesc((tt')); colorbar;
        caxis(sc)
	set(gca,'YDir','normal');
	axis square
        set(gca,'XTick',[1,n])
        set(gca,'XTickLabel',{"0","2\pi"})
        set(gca,'YTick',[1,n])
        set(gca,'YTickLabel',{"0","2\pi"})
	%xlabel('x')
	%ylabel('y')
	%print -depsc thfinal.eps
	
