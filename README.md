# TinkerHap - Read-Based Phasing Algorithm with Integrated Multi-Method Support for Enhanced Accuracy

## Overview
TinkerHap is an accurate read-based phasing tool that integrates multiple methodologies to enhance phasing accuracy. It is designed to efficiently phase genomic data by linking sequencing reads across heterozygous sites, extending haplotype blocks, and incorporating pre-phased data when available.

The full paper detailing the methodology and validation of TinkerHap is available at [bioRxiv](https://doi.org/10.1101/2025.02.16.638517)\
DOI: https://doi.org/10.1101/2025.02.16.638517.

### Key Features:
- **Hybrid Phasing Approach**: Combines read-based phasing with statistical or pedigree-based methods for increased accuracy.
- **Broad Compatibility**: Works with both short-read and long-read sequencing data.
- **High Accuracy**: Achieves superior phasing accuracy for both SNPs and indels.
- **Customizable Outputs**: Generates phased VCF, annotated BAM, and haplotype BED files.

## Requirements
- **Python**: >= 3.6.0
- **Dependencies**:
  - `pysam>=0.17.0`

## Installation
TinkerHap can be installed either by using [pip](https://pip.pypa.io) or directly by cloning the repo and running the python script.
1. Using pip:
```bash
pip install tinkerhap
tinkerhap --help
```

2. Using git: Install pysam, then clone the repository and navigate to the project directory:
```bash
pip install pysam
git clone https://github.com/DZeevi-Lab/TinkerHap.git
cd TinkerHap
python ./tinkerhap.py --help
```

## Usage
TinkerHap is executed via command-line arguments. Below are the available options:

### Required Arguments:
- `-vi` or `--vcf-in`: Input VCF file.
- `-bi` or `--bam-in`: Input BAM file.

### Optional Arguments:
- `-s` or `--sample`: Sample ID to phase.
- `-vs` or `--vcf-scaffold`: Input phased VCF file for scaffold-based phasing.
- `-vo` or `--vcf-out`: Output phased VCF file.
- `-bo` or `--bam-out`: Output BAM file (annotated with phasing information).
- `-ei` or `--bed-in`: Input BED file defining regions to phase.
- `-eo` or `--bed-out`: Output BED file with haplotype blocks.
- `-r` or `--region`: Specific region to phase in the format `RNAME[:STARTPOS[-ENDPOS]]`.
- `-fmq` or `--filter-map-quality`: Minimum mapping quality (default: 20).
- `-md` or `--max-depth`: Maximum read depth to examine (default: 100).
- `-l` or `--log-file`: Log file path for execution logs.
- `-q` or `--quiet`: Quiet mode (suppress console output).

### Example Command:
```bash
python tinkerhap.py -vi input.vcf -bi input.bam -vo phased_output.vcf -bo phased_output.bam -r chr1:100000-200000
```

## Outputs
1. **Phased VCF**: Contains phased variants annotated with phase set (PS) and haplotype information.
2. **Annotated BAM**: Includes phasing information (HP and HT tags).
3. **Split BAM Files**: Separate BAM files for each haplotype.
4. **Haplotype BED File**: Defines haplotype block boundaries for visualization.

## Evaluation
TinkerHap demonstrates high accuracy in phasing both short-read and long-read sequencing data:
- **Short-Read Accuracy**: Up to 96.3% (with hybrid approach).
- **Long-Read Accuracy**: 97.5%.
- **Extended Haplotype Blocks**: Median size of 79,449 base pairs with long-reads.

## Citing

Please cite this paper when using SAMtools for your publications:

> Uri Hartmann, Eran Shaham, Dafna Nathan, Ilana Blech, Danny Zeevi, TinkerHap - A Novel Read-Based Phasing Algorithm with Integrated Multi-Method Support for Enhanced Accuracy, GigaScience, 2025;, giaf138, https://doi.org/10.1093/gigascience/giaf138

## License
TinkerHap is distributed under the MIT License. See `LICENSE` for details.
