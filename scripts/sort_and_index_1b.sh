	#!/bin/bash
	for cnt in {31..48}
	do
	        samtools view -bS RNA150217RH_${cnt}_Aligned.out.sam > RNA150217RH_${cnt}_Aligned.out.bam
	        samtools sort RNA150217RH_${cnt}_Aligned.out.bam  -o RNA150217RH_${cnt}_Aligned_sorted.bam 
	        samtools index RNA150217RH_${cnt}_Aligned_sorted.bam
	done
