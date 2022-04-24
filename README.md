# Folders

- `Datasets` contains all the three datasets used in the experiments.
- `Python` contains the python scripts of EM-HIV.
- `Sample` contains the sample data of training and testing.

# Usage

1. prepare the training and testing datasets by following the format in Sample folder

2. run `python3 EvocleaveV2.py` to extract variable length coevolutionary patterns. Note that EM-HIV will automatically use the variable length coevolutonary patterns extracted from positive set for feature vector construction.

3. run `python3 main.py` to execute EM-HIV. Several parameters have to be predetermined.

   `-f`: the features used t construct feature vectors, possible values of this parameter are 0 (AAI), 1 (CheP), 2 (VLCoP), 3 (AAI+CheP), 4 (CheP+VLCoP), 5 (AAI+VLCoP) and 6 (AAI+CheP+VLCoP);
   `-c1`: the value of C1;
   `-beta`: the value of beta;
   `-i`: input folder;
   `-cv`: optional, PU-HIV will switch to cross validation mode if provided and the value of this parameter is the number of folds in cross validation.
   Hence, a complete command to run the sampl data is `python3 main.py -f 6 -t 30 -c1 8 -beta 2 -i ../Sample`.
   If the paratemter `cv` is provided, the subfolders in the input folder should be named with integers. For example, if the value of `cv` is set as 10, the names of subfolders should be from 1 to 10.

4. check out the results.txt file in the input folder for the prediction results of testing data.

Node: The codes should be compatible with Python 3.7. If you get errors when running the scrips, please try the recommended versions.
