function [x0,y0,x1,y1,x2,y2,x3,y3] = calcular_pos_brazo(thetas, L1, L2, L3)
    th1 = thetas(1); th2 = thetas(2); th3 = thetas(3);
    x0 = 0; y0 = 0;
    x1 = L1*cos(th1); y1 = L1*sin(th1);
    x2 = x1 + L2*cos(th1+th2); y2 = y1 + L2*sin(th1+th2);
    x3 = x2 + L3*cos(th1+th2+th3); y3 = y2 + L3*sin(th1+th2+th3);
end
