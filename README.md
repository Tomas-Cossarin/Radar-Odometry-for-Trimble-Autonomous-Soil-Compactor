# Radar-Odometry-for-Trimble-Autonomous-Soil-Compactor

This project was completed for Trimble in colaboration with a team from the University of Colorado Boulder. The objective was to calculate the change in the autonomous soil compactor's position relative to its starting location, using radar odometry as a very cheap alternative to GNSS.

![compactor](compactor.jpg)

Radar data is received in batches at a variable frequency, typically between 10-18Hz.
These batches are referred to throughout this code as epochs.

At each epoch:
- The data from each radar is rotated to face directly forward
    - The correct rotation for each radar is determined by its known placement on the machine
        and supplemented by the results of a calibration algorithm that fine tunes these adjustments
    
- This adjusted radar data is compiled into two matrices:
    - A matrix consists of [cos(az), sin(az)] for each entry 
    - b matrix consists of radial speeds
    
- A RANSAC implementation from the sklearn library then:
    - Non-deterministically discards outliers in the data
    - Runs a linear regression on the remaining data, i.e., velocity = (A^T A)^-1 (A^T b)
    - The resulting 1x2 velocity vector is in the form [v_x, v_y]
    
- The compactors's location on the x-y plane is updated according
    to this velocity and the length of the epoch in seconds

## Results

Radar odometry was very effective in a straight line.

However, accuracy drops once the vehicle begins turning because the heading is not very accurate.

![ds1x](ds1x.png)
![ds1y](ds1y.png)
![ds2x](ds2x.png)
![ds2y](ds2y.png)
