##Filtering the matrix for only heterozygots sites
vcftools --gzvcf matrix.vcf --extract-FORMAT-info GT | grep "0/1" > filtered_vcf_hete.vcf

##extracting the matrix as text like format
vcftools --vcf filtered_vcf_hete.vcf --extract-FORMAT-info GT --out hete_filtered


##extraction of a specific region based on coordinates
start=5929164
end=5929204
chr=Chr2

vcftools --vcf /groups/nordborg/user/benjamin.jaegle/Documents/001_DATA/101_PAPER_DUPLICATION/for_submission/raw_data/1135g.181k.prior15.ts99.5.BIALLELIC.vcf --chr $chr --from-bp $start --to-bp $end --recode --out $chr.$start.$end

bgzip -c $chr.$start.$end.recode.vcf > $chr.$start.$end.vcf.gz
tabix -f -p vcf $chr.$start.$end.vcf.gz

vcftools --gzvcf $chr.$start.$end.vcf.gz --extract-FORMAT-info GT | grep "0/1"
