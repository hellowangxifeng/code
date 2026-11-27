README for ALFA task 2 M2.2.1 deliverable report September, 2016
Prepared by: Geoffrey Hollinger, Oregon State University, 
geoff.hollinger@oregonstate.edu, 541-737-5906

Sections:
    - Data Description
    - Sample Code
    
************************************************************************
Data Description:

    The data is saved in *.mat files containing the following variables:
        - heading: the heading of the vehicle in degrees
        - t_stamp: the time stamp for the coorisponding other data values
        - x      : the x position of the vehicle
        - y      : the y position of the vehicle
        - z      : the z position of the vehicle
    
    
************************************************************************
Sample Code

    Below is sample code for interacting with the data
    
    Loading:
        - load('P2_stationkeeping_part1_final.mat')
        Using the above command will load the data stored in the corresponding
        *.mat file in the format given above
        
    Calculating Error:
        - [RMSE, ME] = err_calc(actual, desired)
        - [RMSE, ME] = err_calc(x, zeros(length(x)))
        Using the above provided function, err_calc, will give the root mean
        squared error (RMSE) and the mean error (ME) of the given vector
        of measurements to a desired vector of values. The second command can
        be used to get the RMSE and ME of the x position for a data set.
        
        -NOTE: the err_calc() function assumes that the vectors are of the same
               length but will still work is the desired vector is longer than
               the actual vector.
               
    Plotting Position
        - plot3(x,y,z)
        Using the above command will plot the position of the vehicle across
        the trajectory
        
        - plot3(x,y,z,'b'); hold on; plot3(x(1),y(1),z(1),'g-x','LineWidth',3); plot3(x(end),y(end),z(end),'r-o','LineWidth',3)
        Using the above series of commands will plot the position of the vehicle
        across the trajectory and label the start with a green x and the end
        with a red o.    
