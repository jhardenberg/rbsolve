
nt=50000;
for i=1000:1000:nt
    imz(i,129);
    drawnow;
    print('-dpng',sprintf('img%04d.png',i));
end
    