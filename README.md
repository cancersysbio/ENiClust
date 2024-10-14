# ENiClust - Ensemble Integrative Clustering


### How to install?

First clone the git repo:

```bash
git clone https://github.com/cancersysbio/ENiClust.git

cd ENiClust

python setup.py install
```

Once installed, you can type:

```bash
ENiClust.py --help
```


### How to run the classifier with a pre-trained model?


First, preprocess the input data:


```bash
python PreProc.py -g CARIS/eniclust/01_gene_level.txt -s CARIS/eniclust/02_segments.txt -q CARIS/eniclust/03_qc.txt -m CARIS/eniclust/04_mutations.txt -f CARIS/eniclust/05_receptors.txt -y data/general/exome_gene_locations_arminfo.txt -o data/test/output_preproc/caris
```
In this example, `PreProc.py` will preprocess the data set from which the ICs will be predicted, therefore `-r` argument is optional.

Then, run the classifier :


```bash
python ClassifyIC.py -i data/general/training_models/model_44 -d data/test/output_preproc/caris/data_set.txt -o data/test/output_classifyIC/caris
```
In this example, the already trained classifier is stored in `data/general/training_models/model_44`.

Or, all at once :


```bash
python ENiClust.py -g CARIS/eniclust/01_gene_level.txt -s CARIS/eniclust/02_segments.txt -q CARIS/eniclust/03_qc.txt -m CARIS/eniclust/04_mutations.txt -f CARIS/eniclust/05_receptors.txt -y data/general/exome_gene_locations_arminfo.txt -i data/general/training_models/model_44 -o data/test/output_eniclust/caris
```

### Notes
Here, the predicted class “Other” means “IC3/IC7” (see --ic_cluster argument).
