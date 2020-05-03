# PacBio HiFi Read Assembly Snakemake Workflow



## Installation of tools.
The majority of tools used are available via conda environments or as docker images. The following tools were not and require manual installation. The exact commands used are documented. The tools were all installed in a /tools folder on the same level as the pacbio_hifi/ directory. 

### Hifiasm
```
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
```

### Mummer v4.0
```
wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar -xf mummer-4.0.0beta2.tar.gz
rm mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2/
./configure
make
```

### LTR_FINDER_Parallel
```
git clone https://github.com/oushujun/LTR_FINDER_parallel.git
```

### Hicanu
```
git clone https://github.com/marbl/canu.git
cd canu
git checkout hicanu_rc
cd src
make
```
