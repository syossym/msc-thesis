function num_iter = Find_number_iteration(accuracy)
    % Find the number of light round trip to achieve the necessary accuracy
    % in the circulating power
    % simple problem of geometric progression 
    
    global ITM
    global ETM

    % Round trip amplitude gain for the field (= 1 - round trip loss)
    RT_loss = ITM.r*ETM.r;

   % Have to solve RT_loss^num_iter < 0.5*accuracy
    
    num_iter = log(0.5*accuracy)/(log(RT_loss));
`       
end