# GWAS Pipeline for Colorectal Cancer Study

## PLINK QC

## Missing rate per sample, impose 95% call rate
`plink --file PLINK/Smokescreen_Biorealm_p9-10 --mind 0.05 --recode --out results/crc`

## Missing rate per snp, impose 95% call rate
`plink --file results/crc --geno 0.05 --recode --out results/crc`

## Check statistics of the missingness
`plink --file results/crc --missing`

## Filter only MAF > 1%
`plink --file results/crc --maf 0.01 --recode --out results/crc`

## Check statistics of allele frequency
`plink --file results/crc --freq`

## Perform Hardy Weinberg Equilibrium test and report the statistics (p-value < 1e-6)
`plink --file results/crc --hardy midp --hwe 1e-6 midp --recode --out results/crc`

## Check heterozygosity
`plink --file results/crc --het small-sample --out het`

## Create a txt file to save HET information
`echo "FID IID obs_HOM N_SNPs prop_HET" > het.txt`
`awk 'NR>1{print $1,$2,$3,$5,($5-$3)/$5}' het.het >> het.txt`

## Determine 3SD of heterozygosity rates (HR)
`awk 'NR>1{sum+=$5;sq+=$5^2}END{avg=sum/(NR-1);print avg-3*(sqrt(sq/(NR-2)-2*avg*(sum/(NR-2))+(((NR-1)*(avg^2))/(NR-2)))),avg+3*(sqrt(sq/(NR-2)-2*avg*(sum/(NR-2))+(((NR-1)*(avg^2))/(NR-2))))}' het.txt`

## Create a list of samples whose HR values are outside of 3SD range
`awk '$5<=<lower-limit> || $5>= <upper-limit>' het.txt> het.drop`

## Remove the samples 
`plink --file results/crc --remove het.drop --recode --out results/crc`

## Variant LD pruning
`plink --file results/crc --indep-pairwise 50 5 0.5 --out indep`
`plink --file results/crc --extract indep.prune.in --check-sex --out sex2`
`grep PROBLEM sex2.sexcheck > sex.drop`
`grep PROBLEM *.sexcheck`
`plink --file results/crc --remove sex.drop --recode --out results/crc`

## Check relatedness
`plink --file results/crc --chr 1-22 --extract indep.prune.in --genome --out ibd`
`plink --file results/crc --extract indep.prune.in --genome --min 0.2 --out pihat`

## Remove duplicates
`plink --file results/crc --list-duplicate-vars 'ids-only' 'suppress-first' --out results/crc.dupvar`
`plink --file results/crc --exclude results/crc.dupvar --recode --out results/crc`

## Remove indels
`plink --file results/crc --snps-only 'just-acgt' --recode --out results/crc`

## Convert to vcf
`plink --file results/crc --recode vcf --out results/crc`

## Convert to vcf.gz
`bcftools sort results/crc.vcf -Oz -o results/crc.vcf.gz`

## Slice by chromosome
`bcftools index -s results/crc.vcf.gz | cut -f 1 | while read C; do bcftools view -Oz -o results/vcf_per_chr/chr${C}.crc.vcf.gz results/crc.vcf.gz "${C}" ; done`

## Handle the multiallelic sites
`for CHR in {1..22}; do bcftools norm -m -any results/vcf_per_chr/chr${CHR}.crc.vcf.gz -Oz -o results/vcf_per_chr_temp/chr${CHR}.vcf.gz; done`

## Check reference allele mismatches
`bcftools norm --check-ref e -f ~/GRCh37/human_g1k_v37.fasta results/crc.vcf.gz -Ou -o /dev/null`

## Conform-gt for flip strand issues
`for i in {1..22}; do java -jar conform-gt.24May16.cee.jar gt=results/vcf_per_chr/chr$i.crc.vcf.gz ref=~/GRCh37/1kg_per_chr/chr$i.1kg.phase3.v5a.vcf.gz chrom=$i match=POS out=results/vcf_per_chr_conform/crc.chr$i; done`

`bcftools concat -o crc_all.vcf results/vcf_per_chr_conform/chr*.crc.vcf.gz`

`java -jar conform-gt.24May16.cee.jar ref=chr$i.1kg.phase3.v5a.vcf.gz gt=chr20.vcf.gz chrom=20 out=mod.chr20 excludesamples=non.eur.excl`

`bcftools norm --check-ref e -f ~/GRCh37/human_g1k_v37.fasta results/trial.vcf.gz -Ou -o /dev/null`

`for CHR in {1..22}; do bcftools norm -m -any results/vcf_per_chr/chr${CHR}.trial.vcf.gz -Oz -o results/vcf_per_chr_temp/chr${CHR}.split_multiallelic_trial.vcf.gz; done`

`for CHR in {1..22}; do echo ${CHR} chr${CHR}; done >> chr_names.txt`

`for CHR in {1..22}; do bcftools annotate --rename-chrs chr_names.txt results/vcf_per_chr_temp/chr${CHR}trial_split_multiallelic.vcf.gz -Oz -o results/vcf_per_chr_temp/chr${CHR}.trial.vcf.gz; done`

`bcftools concat results/vcf_per_chr_conform/trial.chr{1..22}.vcf.gz -o results/imputed_trial.vcf.gz -Oz`

`for i in {1..22}; do unzip chr_$i.zip && rm chr_$i.zip; done`
