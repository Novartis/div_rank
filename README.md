### Diversity picking code

This code is examplified on the public SSEC library (a library that can be used in screens withe xternal partners at the FASTlab)


The code is run like this:

	python div_rank.py <path to config file>

In this example case:

	python div_rank.py config_curated_25.json


#### Configuration file architecture
This json file defines the input, and output paths and all other runtime parameters. The following paremters are recognized

* **path**: path from which input is read and output is written to
* **rank_file** : a file containg a unique record for each sample, provoding the ranking of each sample, expected to be a csv file. Additionla columns bsides sample ID and rank can be carried through
* **sample_name** : name of the column in the rank_file which contains the sample ID
* **score_column** : column in the rank-file containg the score
* **random_id_column** : The code needs a random number for tie breaking, this can be provided via the rank_file, in which case a column name needs to be given. 
* **target_picking_quorum** : A floting point number in range between 0.0 1nd 1.0. This determines when diversity picking shoud stop. All compounds not yet picked are then added to one last picking round. This meake sense, as tpically a well define dpicking order is not usefull in 
* **outfile** : name of the output file to be written within path
* **class_data**: A dictionary of class data items. Each item represents a classification method. The dictionary key is an abbreviation of the method, whereas the value is a dictionary specifying the necessary setup:
	* **file**: delimited text file with the classification data
    * **sep**: delimiter
    * **min_size**: a list with the minimum class size for each round (for example \[25,5,1\] means that for this method at round 1 class size of >= 25 is required, from round 2 on a size quorum of 5 needs to be met, and from round 3 on singletons are accepted )
    * **sample_name**: name of the column corresponding to the sample name , must match a record in rank_file
    * **class_id**: identifier for the class, does not need to be unique across different classification methods

