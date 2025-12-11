for i in `cat ALDH_samples.txt`
do
    echo -e ">> ${i}"

    bcftools norm -m -both -f ${ref} --threads 15 ${refdir}/${i}.g.vcf.gz -Ou \
    | bcftools sort -Oz -o ${workdir}/${i}.g.norm.vcf.gz

    tabix -@ 15 ${workdir}/${i}.g.norm.vcf.gz
done