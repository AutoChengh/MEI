# Modified-Emergency Index (MEI): Advanced Criticality Metric for Autonomous Driving in Lateral Conflict

## Overview

- This repository provides a codebase designed to evaluate real-time collision risks between two traffic participants based on their motion states and dimensions called **MEI**. 
- First, the Separating Axis Theorem (SAT) is used to determine whether a collision has occurred. For pairs of participants that have not collided, the code calculates two key metrics: **ACT** and **MEI**, which quantify the real-time risk of collision.

## Key Metrics

- **MEI (Minimum Evasive Impulse)**:  
    - Unit: [m/s]  
    - MEI > 0 indicates that a collision is inevitable if the current motion states remain unchanged.  
    - A higher MEI value signifies a higher collision risk.  
    - The goal is to reduce MEI to a value less than 0 to resolve the conflict.

- **ACT (Augmented Collision Time)**:  
    - Unit: [s]  
    - ACT is a 2D extension of the Time-to-Collision (TTC) metric.  
    - ACT > 0, with values closer to 0, indicates higher collision risk.  

This repository also includes an open-source implementation of ACT calculations, specifically for scenarios where MEI > 0 is a prerequisite.

## Library Requirements
This script requires the following Python libraries:
- `numpy`

## Input Data Format

The code processes vehicle motion states using the following variables:

### Variable Descriptions

| **Variable** | **Description** | **Unit** | **Format** |
|------------|----------|----------|----------|
| **x_A** | X-coordinate of ego vehicle (Vehicle A) centroid | meters (m) | float |
| **y_A** | Y-coordinate of ego vehicle (Vehicle A) centroid | meters (m) | float |
| **v_A** | Speed of ego vehicle (Vehicle A) | meters/second (m/s) | float |
| **h_A** | Heading angle of ego vehicle (Vehicle A) | radians (rad) | float, range [-π, π], e.g., 1.57 represents 90° |
| **l_A** | Length of ego vehicle (Vehicle A) | meters (m) | float |
| **w_A** | Width of ego vehicle (Vehicle A) | meters (m) | float |
| **x_B** | X-coordinate of surrounding vehicle (Vehicle B) centroid | meters (m) | float |
| **y_B** | Y-coordinate of surrounding vehicle (Vehicle B) centroid | meters (m) | float |
| **v_B** | Speed of surrounding vehicle (Vehicle B) | meters/second (m/s) | float |
| **h_B** | Heading angle of surrounding vehicle (Vehicle B) | radians (rad) | float, range [-π, π], e.g., 1.57 represents 90° |
| **l_B** | Length of surrounding vehicle (Vehicle B) | meters (m) | float |
| **w_B** | Width of surrounding vehicle (Vehicle B) | meters (m) | float |

## Future Work

- More experiment-related code will be made public in the future. Stay tuned!