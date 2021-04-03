# seqcheck

### Overview

```seqcheck``` is a command line utility to make it easier and faster to check sequencing results.

### The Configuration File

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