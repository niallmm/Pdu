function u0 = initold(xnum, initv)

    h0 = initv(1,:);
    c0 = initv(2,:);

    sizeh0 = size(h0);

    u0(1:2:xnum*2)= h0;

    u0(2:2:xnum*2)= c0;