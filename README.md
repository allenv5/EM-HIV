# Folders

- `Datasets` contains all the three datasets used in the experiments.
- `Evocleave V2.0` contains the variable length coevolutionary patterns from `Sample` folder.
- `Python` contains the python scripts of EM-HIV.
- `Sample` contains the sample data of training and testing.

# Usage

1. prepare the training and testing datasets by following the format in `Sample` folder

2. run `python3 EvocleaveV2.py` to extract variable length coevolutionary patterns. Note that EM-HIV will automatically use the variable length coevolutionary patterns extracted from positive set for feature vector construction. For the convenience of testing, we have uploaded the `Evocleave V2.0` file obtained by running EvocleaveV2.py from the default `Sample` file. If the user needs to test other data, please prepare the data according to the format of the `Sample` file, and execute the command to get the corresponding `Evocleave V2.0` file.

3. run `python3 main.py` to execute EM-HIV. Several parameters have to be predetermined.

   - `-f`: the features used to construct feature vectors, possible values of this parameter are 0 (AAI), 1 (CheP), 2 (VLCoP), 3 (AAI+CheP), 4 (CheP+VLCoP), 5 (AAI+VLCoP) and 6 (AAI+CheP+VLCoP);
   - `-c1`: the value of C1;
   - `-beta`: the value of beta;
   - `-i`: input folder;
   - `-cv`: optional, PU-HIV will switch to cross validation mode if provided and the value of this parameter is the number of folds in cross validation.
   1. In order to facilitate users to run, we have introduced a series of parameter default values. For example, the default value of `f` is 6, the default value of `t` is 30, the default value of `c1` is 8, the default value of `beta` is 2, and the default value of `i` is "../Sample".
     2. Hence, the user only needs to execute the command `python3 main.py` to run the sample data.
     3. If the paratemter `cv` is provided, the subfolders in the input folder should be named with integers. For example, if the value of `cv` is set as 10, the names of subfolders should be from 1 to 10.

4. check out the results.txt file in the input folder for the prediction results of testing data.

Node: The codes should be compatible with Python 3.7. If you get errors when running the scripts, please try the recommended versions.
