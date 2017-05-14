        function im(t,nz)

        filename=sprintf('t%07d.zzz',t);
        tt=load(filename);
        whos
        n=length(tt); n=n/nz;
        tt=reshape(tt,n,nz);
        tt1=tt-ones(n,1)*mean(tt);
         subplot(2,1,1)
        imagesc((tt1')); colorbar;
        set(gca,'YDir','normal');
        xlabel('x');
        ylabel('z');
title(sprintf('t=%d, Re=1e7',t),'FontSize',18)


	print -depsc tvfinal.eps
	
