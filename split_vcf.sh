for f in *.vcf; do
 file=$f

cat  <(bcftools view -h $file | grep ^##)  <(cat header_lines.txt)  <(bcftools view -h $file | grep -v ^##) > header.txt

bcftools reheader --header header.txt $file > test.vcf

file="test.vcf"

 for sample in `bcftools query -l $file`; do
  echo $sample
  echo ${file/.vcf*/.$sample.vcf.gz}
  echo $file
  bcftools view -c1 -Ov -s $sample -o ${file/.vcf*/.$sample.vcf} $file
 done
done
