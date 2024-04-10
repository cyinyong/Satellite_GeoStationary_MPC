function O = Eul2Omat(x)
    O1  = [1 0 0; 0 cos(x(1)) sin(x(1)); 0 -sin(x(1)) cos(x(1))]; 
    O2  = [cos(x(2)) 0 -sin(x(2)); 0 1 0; sin(x(2)) 0 cos(x(2))];
    O3  = [cos(x(3)) sin(x(3)) 0; -sin(x(3)) cos(x(3)) 0; 0 0 1];
    O = O1*O2*O3;   
    
end