# Testing Reference Guide

This file contains the required verification and testing procedures for this project.
All test files and assets shall reside within the test/ folder.

## Unit testing:

Run unit tests for parameter validation and intermediate calculations

- Run all tests: <Insert concrete test commands>
- Run single test: <Insert concrete test commands>

## E2E verification:

1. Run script with projects verification input geometry and parameters. Save the following log from the script for analysis and the pull-request note (TODO: implement as script output in a .log file with the same name as the output file):
- Run data from the script including
	- The runs input parameters 
	- Date and time of run start
	- Duration from start to completion
	- Run characteristics (number of tiles, number of parallel threads per stage of the generation procedure, etc)
	- Maximum memory usage
	
<Insert concrete test commands>	

3. Verification: If no golden sample is defined, provide output file for user verification. If golden sample is defined, check similarity of geometriers by running a geometry check script that inspects if the output geometry is inhibiting the same volume as the golden sample (e.g. subtraction either way should leave near zero volume), ensure that the generated geometries are manifold and non-self-intersecting. 

<Insert concrete test commands>

## Goal oriented performance optimization:

The log should contain sufficient infomration to assist iterative goal oriented performance optimization. 
TODO: Write test and optimization procedure

## Verification Checklist

1. Ensure all linting passes
2. Verify edge cases for any modified boundary logic.
3. [Add testing steps here...]