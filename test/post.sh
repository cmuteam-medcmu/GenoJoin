fname='test_dv'

bgzip ${fname}_h.vcf
tabix ${fname}_h.vcf.gz

gatk --java-options "-Xmx4g" CalculateGenotypePosteriors -V ${fname}_h.vcf.gz -O ${fname}_r.vcf.gz

bcftools view -i 'AC>=2 & AC<=(AN-2)' ${fname}_r.vcf.gz -Ob -o ${fname}_filtered.bcf
bcftools stats ${fname}_filtered.bcf > ${fname}_stats.txt

rm ${fname}_h.vcf.gz* ${fname}_r.vcf.gz