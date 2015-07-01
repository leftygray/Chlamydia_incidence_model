## Instructions for running chlamydia model

From developer: Ewan Cameron <dr.ewan.cameron@gmail.com> 

The main algorithm is contained in abc.run.abc.R which calls other routines to read in the data and run the ABC fitting routine.  

- The fits will eventually stop by themselves once the code judges any further improvements to be outweighed by the time required in computing, but the job may also be killed manually from outside if you lose patience.  
- Nothing is to be gained by running the fitting code for longer than 12 hours, I would expect.  

You will need to make a folder called output to take the output files from abc.run.abc.R. 

- The latest output file number should then be updated in abc.posterior.processing.R and then this script run until it stops (probably another hour).  
- abc.plot.results.R can then be run to make plots.  

Much of the input/output should be able to be deciphered and then updated by someone familiar with the R language and some excel; the code is not so sophisticated that it can take care of all this by itself, so be gentle with it! (and perform sanity checks along the way like testing whether the contents of any given file have actually been read in or not!).  The simulation code may be somewhat indecipherable but I am happy to answer any questions as you go along.  

Hacking this code to do the same job as in our paper for the 2014 and 2015 numbers should be 'easy' enough and worthwhile, ditto for 2016.  But looking forward to 2017 and beyond I would imagine a more advanced solution could (and should) be sought.  One way to do this would be to replace my simulation code with a tailored version of EMOD which supports e.g. the division of the population into key cohorts like MSM or sex worker populations each with their own transmission behaviours and treatment seeking pathways.  I am using EMOD for malaria simulations at the moment and will have a paper out soon describing a calibration solution similar in spirt to that in our chlamydia paper.  Also, I noticed that David is a co-author on Jeff Eaton's recent HIV modelling paper which has applied EMOD successfully in that context.