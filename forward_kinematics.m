function pos = forward_kinematics(thetas, L1, L2, L3)
    th1 = thetas(1);
    th2 = thetas(2);
    th3 = thetas(3);
    
    x = L1*cos(th1) + L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    y = L1*sin(th1) + L2*sin(th1+th2) + L3*sin(th1+th2+th3);
    
    pos = [x; y];
end
