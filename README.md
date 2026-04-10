# LDP4Cycle

# Directory Structure
- cpp/			&emsp;C/C++ codes.
- python/		&emsp;python codes
- data/		&emsp;datasets
- README.md		&emsp;This file.

# Usage

**(1) Install**

Install StatsLib (see cpp/README.md).

Install C/C++ codes as follows.
```
$ cd cpp/
$ make
$ cd ../
```
**(2) (Optional) Download and preprocess Gplus**

Download the [Gplus dataset](https://snap.stanford.edu/data/ego-Gplus.html) and place the dataset in data/Gplus/.

Run the following commands.

```
$ cd python/
$ python3 ReadGPlus.py ../data/Gplus/gplus_combined.txt ../data/Gplus/edges.csv ../data/Gplus/deg.csv
$ cd ../
```

Then the edge file (edges.csv) and degree file (deg.csv) will be output in data/Gplus/.

**(3) (Optional) Download and preprocess IMDB**

Download the [IMDB dataset](https://www.cise.ufl.edu/research/sparse/matrices/Pajek/IMDB.html) and place the dataset in data/IMDB/.

Run the following commands.

```
$ cd python/
$ python3 ReadIMDB.py ../data/IMDB/IMDB.mtx ../data/IMDB/edges.csv ../data/IMDB/deg.csv
$ cd ../
```
**(4) Run subgraph counting algorithms**

Run the following commands ([Dataset] is "Gplus" or "IMDB" or any datasets listed in data/).

```
$ cd cpp/
$ chmod +x run
$ ./run [Dataset (Gplus/IMDB/etc.)]
$ cd ../
```
For more details of parameter settings, see Usage of SubgraphShuffle.

**(5) (Optional) Visualize the data**

Run the following commands ([Dataset] is "Gplus" or "IMDB" or any datasets listed in data/).

```
$ cd data/
$ python3 plot_results.py
$ cd ../
```


# Execution Environment
We used CentOS 7.5 with gcc 4.8.5.
