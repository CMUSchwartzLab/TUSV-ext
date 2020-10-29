for patient in `aws s3 ls s3://tusv-xuecong-191122-east/data_full/input/`;
	do patient=${patient%/*}
	echo ${patient}
	if [ "${patient}" == "01_TCGA-B6-A0RT" ] || [ "${patient}" == "PRE" ];
	then 
		continue;
	fi	
	mkdir ../data/output/${patient};
	aws s3 cp s3://tusv-xuecong-191122-east/data_full/input/${patient}  ../data/input/${patient} --recursive;
	rm ../data/input/${patient}/*.maf
	cat ../data/input/${patient}/sample1.vcf ../data/input/${patient}/sampled_snvs.vcf > ../data/input/${patient}/sample1_snv.vcf
	rm ../data/input/${patient}/sample1.vcf
	rm ../data/input/${patient}/sampled_snvs.vcf
	rm ../data/input/${patient}/snvs.vcf
	python tusv_1130.py -i ../data/input/${patient} -o ../data/output/${patient} -n 4 -c 10 -t 4 -r 2 -p 1 -m 7200 -s 0 -b > ../out/${patient};
	wait;
	rm -r ../data/input/${patient};
	done
