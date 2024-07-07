# CBS-ML
Machine learning-guided conflict selections for Conflict-Based Search.


## Setup

Follow instructions in [CBSH2](https://github.com/Jiaoyang-Li/CBSH2-RTC) to compile the c++ code.

A set of instances for the lak503d map are in the instance folder.


## Data Collection and Training

Below is an example to collect data for the lak503d map for the instance ```instances/lak503dmap-100agents-0.agents``` with 2 hours runtime. The option ```-u 1``` specifies data collection mode.

```
./CBSH2 -m instances/lak503d.map -o collect_log.csv -t 7200 -s 1 -h WDG -a instances/lak503dmap-100agents-1.agents -u 1 --feature feature/feature_Lak503d_0.txt
```

After data collection, process the data with the following python script for it to be compatible with the SVM rank software.

```
python dataProcessor.py
```

After that, one can use the [SVM rank software](https://www.cs.cornell.edu/people/tj/svm_light/svm_rank.html) to train the model using the produced data files. 

Once finished training, there will be a trained SVM rank model the main folder. Open the file and extract the weights into a new weight file ```weight_file``` with the content being rows of pairs ```featureID weights```.


## Testing

To run test on the lak503d map with 1 hour runtime on instance ```instances/lak503dmap-100agents-2.agents```, one can use the command below. 
The option ```-u 4``` specifies test mode.

```
./CBSH2 -m city/Paris/lak503d.map -o run_test_ML.csv -t 3600 -s 1 -h WDG -a instances/lak503dmap-100agents-2.agents -k 100 -u 4  --model weight_file
```



# Citation

Please cite our paper if you use this code in your work.

```
@inproceedings{huang2021conflict,
  title={Learning to resolve conflicts for multi-agent path finding with conflict-based search},
  author={Huang, Taoan and Dilkina, Bistra and Koenig, Sven},
  booktitle={AAAI Conference on Artificial Intelligence},
  year={2021},
  pages={11246--11253}
}
```
