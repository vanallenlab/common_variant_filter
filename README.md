# common_variant_filter
A simple filter to annotate and filter somatic variants based on their occurence in [ExAC](http://exac.broadinstitute.org/). This is similar to the filter applied in [VCF2MAF's common variant filter](https://github.com/mskcc/vcf2maf/blob/master/docs/vep_maf_readme.txt), as used by [AACR Project GENIE](http://cancerdiscovery.aacrjournals.org/content/7/8/818). A variant is filtered if at least 10 alleles containing that variant are present across any subpopulation in ExAC, unless it appears at a [known somatic site](https://github.com/mskcc/vcf2maf/blob/v1.6.12/data/known_somatic_sites.bed). 

This filter was originally implemented by [Cyriac Kandoth](https://github.com/ckandoth) for [VCF2MAF](https://github.com/mskcc/vcf2maf). 

## Getting common_variant_filter
This codebase is available for download through this Github repository, the [v1.0.0 release](https://github.com/vanallenlab/phial/releases), and [Dockerhub](https://hub.docker.com/r/vanallenlab/phial/).

To download via Github
```
git clone https://github.com/vanallenlab/common_variant_filter
```

To download via Dockerhub
```
docker pull vanallenlab/common_variant_filter:0.1.0
```

## Using common_variant_filter
common_variant_filter accepts the following arguments
- `maf`: A path to a MAF file containing somatic variants for annotation and filtration
- `exac`: A path to ExAC, formatted with the required columns. See [build_exac.sh]() in the `exac` folder. 
- `whitelist`: (Optional) Default `False`. A path to a whitelist of known somatic events in a bed file format. If not passed, no rescuing of somatic sites will occur.
- `min_exac_ac`: (Optional) Default `10`. Minimum allele count across any population to be considered a germline variant. 
- `filter_syn`: (Optional) Default `False`. Filters to the variant classifications of `Missense_Mutation`, `Nonsense_Mutation`, and `Splice_Site`. 
- `min_depth`: (Optional) Default `0`. Filters based on a specified depth of coverage. 

An example run command would be the following,
`python common_variant_filter.py --exac datasources/exac.lite-pass1.4-r1.txt --maf /path/to/maf/sample.maf --whitelist datasources/known_somatic_sites.bed` 

## References
1. [Lek M, Karczewski KJ, Minikel EV, et al. Analysis of protein-coding genetic variation in 60,706 humans. Nature. 2016;536(7616):285-91.](https://www.nature.com/articles/nature19057)
2. [AACR Project GENIE: Powering Precision Medicine through an International Consortium. Cancer Discov. 2017;7(8):818-831.](http://cancerdiscovery.aacrjournals.org/content/7/8/818)
