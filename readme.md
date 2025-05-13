# README

## Overview

This repository provides a codebase designed to evaluate real-time collision risks between two traffic participants based on their motion states and dimensions. The core logic utilizes the Separating Axis Theorem (SAT) to determine whether a collision has already occurred. For pairs of participants that have not collided, the code calculates two key metrics: ACT and MEI, which quantify the real-time risk of collision.

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
