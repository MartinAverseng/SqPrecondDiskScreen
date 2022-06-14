clear all
close all
delete(gcp('nocreate'));

levels = 1:4;

for b1 = 0:1
    for b2 = 0:1
        for b3 = 0:1
            job(levels,1 - b1,1 - b2,1 - b3)
        end
    end
end
