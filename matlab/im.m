	function im(t)

	filename=sprintf('t%07d.dat',t);
        tt=load(filename); 
	n=length(tt); n=sqrt(n);
	tt=reshape(tt,n,n);
	tt=tt-mean(mean(tt));
	imagesc((tt')); colorbar;
	set(gca,'YDir','normal');
	axis square
	%xlabel('x')
	%ylabel('y')
	print -depsc thfinal.eps
	
