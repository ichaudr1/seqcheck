# seqcheck

### Overview

```seqcheck``` is a command line utility to make it easier and faster to check sequencing results.

### Installation

#### You must have Python 3.X.X installed for this to work. 

Installation has only been automates for MacOS. See below for general requirements for Windows and Linux (or MacOS if your prefer to set it up manually or have problems with the automatic install).

Run ```sh setup```. If no errors are present, you should be able to run the utility by calling ```seqcheck``` in terminal in any directory. Type ```seqcheck -h``` for information on usage or see below for more information. 

### Usage

The utility is called by typing ```seqcheck```. 

To check sequencing results, run ```seqcheck --path=[path to directory with results]```. Make sure the construct library path is defined in the configuration file (see below). For each ```.seq``` file in the results directory, you will be asked to enter a construct to compare the sequence against. If the construct you enter does not exist exactly as entered, you will be presented with a list of constructs that are similar to the entered construct name.

Note: If your sequencing results do not have ```.seq``` files in the directory, note an Issue on this repository or email (ichaudr1@umbc.edu). 

To add a new construct(s) to the construct library, use ```seqcheck --add_construct```. This will ask for a path to a csv file with the new constructs in this format: ```construct_name,sequence```. These will be added to the reference construct library that is listed in the configuration file (see below). 

### The Configuration File

This file can be updated by running ```seqcheck --configure``` or it is stored at ```$HOME/bin/seqcheck/config.json``` if installed with the automated script.

```json
{
	"alignment_parameters":[{

		"match_score": 1,
		"mismatch_score": -1,
		"open_gap_score": -1,
		"extend_gap_score": -0.5,
		"target_end_gap_score": 0,
		"query_end_gap_score": 0,
		"max_alignments": 1
	}],

	"construct_lib":"/Users/ichaudr/bin/seqcheck/construct_lib.csv",

	"construct_search_threshold":10
}

```

The ```alignment parameters``` will work as it, but can be changed if needed. More info about these penalties here: [More info](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#Scoring_systems)

The ```construct_lib``` is the path the .csv file holding a list of your constructs and their sequences.

```construct_search_threshold``` is a parameter that defines the stringency of the search function when looking for a construct in the construct library. A lower value is more stringent and a higher value is less stringent.

### Dependencies and Manual Install

Download this repository to a location where you can leave it permanently. 

All dependencies are listed in the ```requirements.txt``` file. You must have Python 3.X.X and a python environment set up with those dependencies. An ```conda``` or ```virtualenv``` will work. The specific instructions will vary based on platform. Once an evironemnt is set up, run ```pip3 install -r requirements.txt```. 

Add the current directory to your path so that the ```seqcheck``` can be ran from terminal outside of the current directory. This will also vary based on platform, but for MacOS, run ```export PATH=$PATH:[current_directory]```. 