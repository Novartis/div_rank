### Diversity picking algorithm

This code allows is for diversity picking across multiple different, and potentially overlapping chemcial compound classes, while at te same time optimizing a property score. 
This algorithm has been used in the re-design of the Novartis screening deck as described in: 

Evolution of Novartis’ Small Molecule Screening Deck Design, Schuffenhauer, A. Schneider N. et al. *J. Med.Chem.* **2020**, https://dx.doi.org/10.1021/acs.jmedchem.0c01332 


Selection algorithm
1.	Initialization:

	* Initialize the selection round with one
	* Determine for each class the initial round, that is the selection round from which onwards it is available for selection based on its membership count and the count limit set for its classification type.  
	* Mark all compounds as not selected.
	* Assign each compound a random number for tie breaking
2.	For each class available for selection at selection round, determine the property rank score of the best ranking compound that has not yet been selected to determine the selection threshold for each class at this round. Skip classes without any compound left for selection.
3.	For each class, count the number of compounds already selected with a property rank score better or equal to selection threshold in previous rounds. Skip the subsequent steps for each class, where this count is ≥ selection round-initial-round+1. This represents the expected number of compounds selected for this class at selection round in case it had no overlap with overlap classes.  
4.	For each remaining class determine the list of selection candidates consisting of all compounds having a  property rank score better or equal to selection threshold 
5.	If there is only one candidate for this class, this is the compound to be selected for the class at this round. If there are multiple selection candidates, select the compound for this class as follows:

	* Choose the compound having the maximal number of class assignments in different class types. 
	* Use the random number assigned at initialization to break remaining ties
6.	Remove the replicates from the list of compounds selected at 5 and mark the compounds as selected in selection round. 
7.	Increment the selection round counter by one 
8.	Repeat steps 2-7 until all compounds are selected. Optionally stop early after a predefined number of selection rounds or if given number of compounds have been selected.   

This code is examplified on the public SSEC library (a library that can be used in screens withe external collaboration partners at the Novartis FASTlab)

### Executing the code
This code is meant to be executed with python 3.6 and pandas 0.25. It has *NOT* been tested under pandas 1.0. 

#### invocation via command line

The code is run like this:

	python div_rank.py <path to config file>

In this example case:

	python div_rank.py config_curated_25.json


#### Configuration file architecture
This json file defines the input, and output paths and all other runtime parameters. The following parameters are recognized:

- **path**: path from which input is read and output is written to. 
- **rank_file** : a file containg a unique record for each sample, provoding the ranking of each sample, expected to be a csv file. Additionla columns bsides sample ID and rank can be carried through
- **sample_name** : name of the column in the rank_file which contains the sample ID
- **score_column** : column in the rank-file containg the score
- **random_id_column** : The code needs a random number for tie breaking, this can be provided via the rank_file, in which case a column name needs to be given. 
- **target_picking_quorum** : A floting point number in range between 0.0 1nd 1.0. This determines when diversity picking shoud stop. All compounds not yet picked are then added to one last picking round. This makes sense, as typically a well define dpicking order is not usefull in towards the end of the picking cycle, as there is litzle benefit in screening for example 90% of the deck versus the full deck
- **outfile** : name of the output file to be written within path
- **class_data**: A dictionary of class data items. Each item represents a classification method. The dictionary key is an abbreviation of the method, whereas the value is a dictionary specifying the necessary setup:
	- **file**: delimited text file with the classification data
    - **sep**: delimiter
    - **min_size**: a list with the minimum class size for each round (for example \[25,5,1\] means that for this method at round 1 class size of >= 25 is required, from round 2 on a size quorum of 5 needs to be met, and from round 3 on singletons are accepted )
    - **sample_name**: name of the column corresponding to the sample name , must match a record in rank_file
    - **class_id**: identifier for the class, does not need to be unique across different classification methods

#### Example data for the SSEC library
The example files provided are covering those screening compounds to which Novartis can give access within collobrations through the FASTLab (**F**acilitaded **A**ccess to **S**creening **T**echnology) program:
* **mol_rank_data.csv**: The ranked compounds which serve as input to the diversity ranking. This contains the finer grained *score* column controlling the preference for picking, the discrete *column_category* giving the property column category. 
* **config_curated_25.json**: The configuration file controlling the picking process. In ordere to run this example you will have to adjust the *path* attribute in this file and let point to the folder where you have cloned this repositoty.
* **cluster_data/\***: The files contained the clsutering information. In the target annotation classes the target name has been deidentified.
* **div_rank_ssec_results_25.csv**: the resulting diversity picking containg the columns present also in the input plus the *picking round* column indicating at which picking round the compound was selected 
