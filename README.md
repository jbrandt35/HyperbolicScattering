# HyperbolicSims


Directed to Josh:

Here's the code I've been using for the last couple months to generate large sets of simulations with chosen parameters. Most of the original infrastructure you wrote relating to submitting large numbers of jobs on PACE is still there. The only thing you would have to change is the names of the directories. To start a simulation set, just run build.sh in the parent directory. The things I've added:

- Change parameters in Setup.py. I've added more options here and it now creates a json file for each set to keep track of the parameters. This can get a bit confusing because of all the options. For example, if you choose randomize_mass = True, then the inner_mass and perturber_mass values will be irrelvent, those are for specific selections. 
- When build.sh is run, it moves all data in the `runs/current_sims` directory to a new directory `runs/{date and time}` to prevent accidental overwriting. build.sh then creates a new directory called `current_sims` where the new data will be saved.
- I just changed how the minimum distance in each run is being saved, so keep an eye out if that causes issues, it might have bugs.

Anything else unclear or not working, let me know, I have things fairly optimized at this point to just running a single command to start everything up.
