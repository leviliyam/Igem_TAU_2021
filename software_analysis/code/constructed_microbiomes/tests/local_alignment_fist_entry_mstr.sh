# !/bin/sh
  cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests/test_local_align_KDVY_first_entry/
chmod 777 ./*
qsub -q TullerNano -r y  -e 0_job_error.txt -o 0_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 0_start_0_end_250_job.sh
qsub -q TullerNano -r y  -e 1_job_error.txt -o 1_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 10_start_2500_end_2750_job.sh
qsub -q TullerNano -r y  -e 2_job_error.txt -o 2_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 11_start_2750_end_3000_job.sh
qsub -q TullerNano -r y  -e 3_job_error.txt -o 3_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 12_start_3000_end_3250_job.sh
qsub -q TullerNano -r y  -e 4_job_error.txt -o 4_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 13_start_3250_end_3500_job.sh
qsub -q TullerNano -r y  -e 5_job_error.txt -o 5_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 14_start_3500_end_3750_job.sh
qsub -q TullerNano -r y  -e 6_job_error.txt -o 6_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 15_start_3750_end_4000_job.sh
qsub -q TullerNano -r y  -e 7_job_error.txt -o 7_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 16_start_4000_end_4250_job.sh
qsub -q TullerNano -r y  -e 8_job_error.txt -o 8_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 17_start_4250_end_4500_job.sh
qsub -q TullerNano -r y  -e 9_job_error.txt -o 9_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 18_start_4500_end_4750_job.sh
qsub -q TullerNano -r y  -e 10_job_error.txt -o 10_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 19_start_4750_end_5000_job.sh
qsub -q TullerNano -r y  -e 11_job_error.txt -o 11_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 1_start_250_end_500_job.sh
qsub -q TullerNano -r y  -e 12_job_error.txt -o 12_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 20_start_5000_end_5250_job.sh
qsub -q TullerNano -r y  -e 13_job_error.txt -o 13_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 21_start_5250_end_5500_job.sh
qsub -q TullerNano -r y  -e 14_job_error.txt -o 14_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 22_start_5500_end_5750_job.sh
qsub -q TullerNano -r y  -e 15_job_error.txt -o 15_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 23_start_5750_end_6000_job.sh
qsub -q TullerNano -r y  -e 16_job_error.txt -o 16_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 24_start_6000_end_6250_job.sh
qsub -q TullerNano -r y  -e 17_job_error.txt -o 17_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 25_start_6250_end_6500_job.sh
qsub -q TullerNano -r y  -e 18_job_error.txt -o 18_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 26_start_6500_end_6750_job.sh
qsub -q TullerNano -r y  -e 19_job_error.txt -o 19_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 27_start_6750_end_7000_job.sh
qsub -q TullerNano -r y  -e 20_job_error.txt -o 20_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 28_start_7000_end_7250_job.sh
qsub -q TullerNano -r y  -e 21_job_error.txt -o 21_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 29_start_7250_end_7500_job.sh
qsub -q TullerNano -r y  -e 22_job_error.txt -o 22_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 2_start_500_end_750_job.sh
qsub -q TullerNano -r y  -e 23_job_error.txt -o 23_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 30_start_7500_end_7750_job.sh
qsub -q TullerNano -r y  -e 24_job_error.txt -o 24_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 31_start_7750_end_8000_job.sh
qsub -q TullerNano -r y  -e 25_job_error.txt -o 25_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 32_start_8000_end_8250_job.sh
qsub -q TullerNano -r y  -e 26_job_error.txt -o 26_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 33_start_8250_end_8500_job.sh
qsub -q TullerNano -r y  -e 27_job_error.txt -o 27_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 34_start_8500_end_8750_job.sh
qsub -q TullerNano -r y  -e 28_job_error.txt -o 28_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 35_start_8750_end_9000_job.sh
qsub -q TullerNano -r y  -e 29_job_error.txt -o 29_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 36_start_9000_end_9250_job.sh
qsub -q TullerNano -r y  -e 30_job_error.txt -o 30_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 37_start_9250_end_9500_job.sh
qsub -q TullerNano -r y  -e 31_job_error.txt -o 31_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 38_start_9500_end_9750_job.sh
qsub -q TullerNano -r y  -e 32_job_error.txt -o 32_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 39_start_9750_end_10000_job.sh
qsub -q TullerNano -r y  -e 33_job_error.txt -o 33_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 3_start_750_end_1000_job.sh
qsub -q TullerNano -r y  -e 34_job_error.txt -o 34_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 40_start_10000_end_10250_job.sh
qsub -q TullerNano -r y  -e 35_job_error.txt -o 35_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 41_start_10250_end_10500_job.sh
qsub -q TullerNano -r y  -e 36_job_error.txt -o 36_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 42_start_10500_end_10750_job.sh
qsub -q TullerNano -r y  -e 37_job_error.txt -o 37_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 43_start_10750_end_11000_job.sh
qsub -q TullerNano -r y  -e 38_job_error.txt -o 38_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 44_start_11000_end_11250_job.sh
qsub -q TullerNano -r y  -e 39_job_error.txt -o 39_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 45_start_11250_end_11500_job.sh
qsub -q TullerNano -r y  -e 40_job_error.txt -o 40_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 46_start_11500_end_11750_job.sh
qsub -q TullerNano -r y  -e 41_job_error.txt -o 41_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 47_start_11750_end_12000_job.sh
qsub -q TullerNano -r y  -e 42_job_error.txt -o 42_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 48_start_12000_end_12250_job.sh
qsub -q TullerNano -r y  -e 43_job_error.txt -o 43_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 49_start_12250_end_12396_job.sh
qsub -q TullerNano -r y  -e 44_job_error.txt -o 44_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 4_start_1000_end_1250_job.sh
qsub -q TullerNano -r y  -e 45_job_error.txt -o 45_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 5_start_1250_end_1500_job.sh
qsub -q TullerNano -r y  -e 46_job_error.txt -o 46_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 6_start_1500_end_1750_job.sh
qsub -q TullerNano -r y  -e 47_job_error.txt -o 47_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 7_start_1750_end_2000_job.sh
qsub -q TullerNano -r y  -e 48_job_error.txt -o 48_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 8_start_2000_end_2250_job.sh
qsub -q TullerNano -r y  -e 49_job_error.txt -o 49_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 9_start_2250_end_2500_job.sh
