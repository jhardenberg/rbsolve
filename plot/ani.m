function ani(fn, t1, t2, step, speed)
    ns=(t2-t1)/step+1;

    if ~exist('speed')
       speed=20;
    end

    f = figure('visible','off','Position',[0, 0, 1408 384]);
    v = VideoWriter(fn,'MPEG-4');
    v.FrameRate = speed;
    open(v)
    imtz(t1)
    set(gca,'Position',[0.035 0.1150 0.91 0.8150])
    for i=t1:step:t2
        fprintf('%i\n', i)
        imtz(i);
        frm = getframe(gcf);
        writeVideo(v, frm)
    end
    close(v)
    close(f)
end



