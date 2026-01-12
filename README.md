# Plant Cas9 CRISPRa Guide Design Pipeline

Specific challenges solved:
- Run extensive off-target analysis using [Sassy](https://github.com/RagnarGrootKoerkamp/sassy)
- Score guide efficiency using [RS3](https://github.com/gpp-rnd/rs3/blob/master/README.md)
- Prioritize guides based on [CRISPR-Act3.0](https://www.nature.com/articles/s41477-021-00953-7) paper and RS3


## Installation

Clone this repo
```{bash}
$ git clone https://github.com/samuelmf1/plant-cas9-crispra.git
$ cd plant-cas9-crispra
```

Create Conda environment
```{bash}
$ conda env update -f environment.yaml --prune
```

If you run into issues with Sassy, install it manually (you will need [Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html))
```{bash}
$ git clone https://github.com/RagnarGrootKoerkamp/sassy.git
$ cd sassy
$ RUSTFLAGS="-C target-cpu=native" cargo install --path . --force
```

## Usage

```{bash}
$ nextflow run main.nf \
    -profile slurm \
    --genome data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    --annotation data/annotation.csv

Nextflow 25.10.2 is available - Please consider updating your version to it

 N E X T F L O W   ~  version 25.10.0

Launching `main.nf` [fervent_kay] DSL2 - revision: d6cd74077d

executor >  slurm (204)
[0d/82804a] PROMOTER_GUIDES (1) | 1 of 1 ✔
[66/8986d2] SPLIT (1)           | 1 of 1 ✔
[b6/9a539f] OFFTARGETS (25)     | 100 of 100 ✔
[94/1e8782] RANK_GUIDES (100)   | 100 of 100 ✔
[15/fd87c9] MERGE               | 1 of 1 ✔
[a8/70d25e] SELECT_GUIDES       | 1 of 1 ✔
Completed at: 12-Jan-2026 20:30:48
Duration    : 10m 38s
CPU hours   : 9.9
Succeeded   : 204
```