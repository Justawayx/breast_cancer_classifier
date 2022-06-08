#PBS -N DESeq2
#PBS -o /mnt/BioAdHoc/Groups/vd-vijay/vfajardo/BISB/BENG-203/course_project/jobs_scripts/feature_engineering/DEA_deseq2_DiseaseStatus_All.out.txt
#PBS -e /mnt/BioAdHoc/Groups/vd-vijay/vfajardo/BISB/BENG-203/course_project/jobs_scripts/feature_engineering/DEA_deseq2_DiseaseStatus_All.err.txt
#PBS -m abe
#PBS -M vfajardo@lji.org
#PBS -q default
#PBS -l nodes=1:ppn=2
#PBS -l mem=50gb
#PBS -l walltime=24:00:00

echo -e "\n\n############    --   Differential Expression Analysis   --    ############"
echo -e "############    -----------  Based on DESeq2   -----------    ############\n\n"

Rscript /home/vfajardo/scripts/differential_analysis/deseq2_based/DESeq2_based_DEA.1.0.R --CountData /mnt/hpcscratch/vfajardo/BISB/BENG-203/course_project/data/MergedRawData.csv --MetaData /mnt/hpcscratch/vfajardo/BISB/BENG-203/course_project/data/MergedMetadata.csv --ReportsPath /mnt/BioAdHoc/Groups/vd-vijay/vfajardo/BISB/BENG-203/course_project/feature_engineering/deseq2_DiseaseStatus_All --Design "~ disease.status"  --LFCShrink TRUE --PThold 1.1 --LFCThold 0 --TPMData /mnt/hpcscratch/vfajardo/BISB/BENG-203/course_project/data/MergedTPMData.csv

echo -e "Job completed!\nCheck for errors if any."
