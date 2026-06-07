function thetas_final = inverse_kinematics_solver(target_pos, current_thetas, L1, L2, L3)
    thetas = current_thetas;
    tol = 0.01;      
    max_iter = 100;  
    
    for i = 1:max_iter
        current_pos = forward_kinematics(thetas, L1, L2, L3);
        err = target_pos - current_pos;
        
        if norm(err) < tol
            break; 
        end
        
        J = jacobian_matrix(thetas, L1, L2, L3);
        thetas = thetas + pinv(J) * err;
    end
    
    thetas_final = thetas;
end
