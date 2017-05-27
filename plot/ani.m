
nt=50000;
for i=0:500:nt
    imtz(i,192,[]);
    drawnow;
    print('-dpng',sprintf('img%05d.png',i));
end
    