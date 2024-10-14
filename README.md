# ENiClust - Ensemble Integrative Clustering
> IC-subtypes prediction from whole-exome sequencing (WES) or genome genome sequencing (WGS)

## Installation

To install ENiClust, follow these steps:

1. Clone the Git repository:

```bash
   git clone https://github.com/cancersysbio/ENiClust.git
```

2. Change to the ENiClust directory:

```bash
cd ENiClust
```
3. Run the setup script to install:
```bash
python setup.py install
```

Once installed, you can view the available commands and options by typing:

```bash
ENiClust.py --help
```



## Running the Classifier with a Pre-trained Model

### Step 1: Preprocess the Input Data

Use `PreProc.py` to preprocess your input data files:


```bash
python PreProc.py -g CARIS/eniclust/01_gene_level.txt \
                  -s CARIS/eniclust/02_segments.txt \
                  -q CARIS/eniclust/03_qc.txt \
                  -m CARIS/eniclust/04_mutations.txt \
                  -f CARIS/eniclust/05_receptors.txt \
                  -y data/general/exome_gene_locations_arminfo.txt \
                  -o data/test/output_preproc/caris
```


In this example:

 `PreProc.py` preprocesses the input data for clustering.
    The `-r` argument is optional when predicting ICs (Integrative Clusters).


### Step 2: Run the Classifier

Run the classifier with the pre-trained model:

```bash
python ClassifyIC.py -i data/general/training_models/model_44 \
                     -d data/test/output_preproc/caris/data_set.txt \
                     -o data/test/output_classifyIC/caris
```

In this example:
> The pre-trained model is stored in `data/general/training_models/model_44`.

### Optional: Run All Steps at Once

Alternatively, you can run the preprocessing and classification steps together:

```bash
python ENiClust.py -g CARIS/eniclust/01_gene_level.txt \
                   -s CARIS/eniclust/02_segments.txt \
                   -q CARIS/eniclust/03_qc.txt \
                   -m CARIS/eniclust/04_mutations.txt \
                   -f CARIS/eniclust/05_receptors.txt \
                   -y data/general/exome_gene_locations_arminfo.txt \
                   -i data/general/training_models/model_44 \
                   -o data/test/output_eniclust/caris
```


### Notes
Here, the predicted class “Other” means “IC3/IC7” (see --ic_cluster argument).
