# What is this and how do I use it

This sub-tree is written by Darach Miller, it takes a zipped up data archive
(in pieces, [in this OSF repo](https://osf.io/7yt59/)). Download that.
Then you can stick
those back together and unzip it to make the `data` directory by running

    make unfold_archive

Then on a modern system running `Singularity`, you can do `make nextflow` to
`curl` down the `nextflow` script, and that'll handle downloading containers
from Singularity hub and running them. Those containers are 'frozen' as of
now. It should run everything, and take up ~15GB of disk space in the work
directory.

Consult in `reports` if you'd like examples of the run times, resource usage,
a DAG of the workflow. It's a bit messy because I'm using `into` and other
operators not as they are intended to be used...

This is a great launching part if you want to dust off the accumulation curve
analysis, the homodimer analysis, or if you just want a few sqlite databases
with the counts and fitnesses in there.

# Running

Change the `scripts/run_pipeline.nfconfig` file to reflect your system
capabilities ! This was re-run finally before submitting the paper on a 
Dell R710 with 16 cores and 72GB RAM, so adjust for your computer !
Look at the 'archtypal' report in the `reports` folder for more resource
requirement info.

# Design

The automation is all controlled at the top level by `make`. This gets and
launches `nextflow` pipelines as appropriate. There's a few other `Makefile`
rules to do other stuff, clean up and archiving data (or unarchiving).

If you want to do custom stuff, I recommend you run it once, then see where the
`nextflow` has made temporary directories in `work`, then you can change in 
there and you'll have sym links to all the raw materials for that step.

# Structure

Makefile:

Use this to actually do things. Use those rules. 

data:

Is raw and or processed data that would be distributable so someone may re-run 
it

tmp:

Is intermediate files

scripts:

Are scripts and code 

reports:

Reports from nextflow, including a DAG ! Neato !
