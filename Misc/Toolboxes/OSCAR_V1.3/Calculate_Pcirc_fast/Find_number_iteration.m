function num_iter = Find_number_iteration(accuracy)
    % Find the number of light round trip to achieve the necessary accuracy
    % (in power)
    % in the results
    % simple problem of geometric progression 
    % accuracy = 0.01 for 1% accuracy
    
    global ITM
    global ETM

    % Round trip amplitude gain for the field (= 1 - round trip loss)
    RT_loss = ITM.r*ETM.r;

   % Have to solve RT_loss^num_iter < 0.5*accuracy
    
    num_iter = floor(log(0.5*accuracy)/(log(RT_loss)));
       
end