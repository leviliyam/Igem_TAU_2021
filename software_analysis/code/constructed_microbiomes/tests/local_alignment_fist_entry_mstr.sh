# !/bin/sh
 cd / tamir1 / liyamlevi / projects / communique / Igem_TAU_2021 / software_analysis / code / raw_data_analysis
chmod 777 ./*
qsub -q TullerNano -r y  -e 0_job_error.txt -o 0_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 0_start=0_end=250.fasta_exec.py
qsub -q TullerNano -r y  -e 1_job_error.txt -o 1_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 10_start=2500_end=2750.fasta_exec.py
qsub -q TullerNano -r y  -e 2_job_error.txt -o 2_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 11_start=2750_end=3000.fasta_exec.py
qsub -q TullerNano -r y  -e 3_job_error.txt -o 3_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 12_start=3000_end=3250.fasta_exec.py
qsub -q TullerNano -r y  -e 4_job_error.txt -o 4_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 13_start=3250_end=3500.fasta_exec.py
qsub -q TullerNano -r y  -e 5_job_error.txt -o 5_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 14_start=3500_end=3750.fasta_exec.py
qsub -q TullerNano -r y  -e 6_job_error.txt -o 6_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 15_start=3750_end=4000.fasta_exec.py
qsub -q TullerNano -r y  -e 7_job_error.txt -o 7_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 16_start=4000_end=4250.fasta_exec.py
qsub -q TullerNano -r y  -e 8_job_error.txt -o 8_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 17_start=4250_end=4500.fasta_exec.py
qsub -q TullerNano -r y  -e 9_job_error.txt -o 9_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 18_start=4500_end=4750.fasta_exec.py
qsub -q TullerNano -r y  -e 10_job_error.txt -o 10_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 19_start=4750_end=5000.fasta_exec.py
qsub -q TullerNano -r y  -e 11_job_error.txt -o 11_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 1_start=250_end=500.fasta_exec.py
qsub -q TullerNano -r y  -e 12_job_error.txt -o 12_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 20_start=5000_end=5250.fasta_exec.py
qsub -q TullerNano -r y  -e 13_job_error.txt -o 13_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 21_start=5250_end=5500.fasta_exec.py
qsub -q TullerNano -r y  -e 14_job_error.txt -o 14_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 22_start=5500_end=5750.fasta_exec.py
qsub -q TullerNano -r y  -e 15_job_error.txt -o 15_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 23_start=5750_end=6000.fasta_exec.py
qsub -q TullerNano -r y  -e 16_job_error.txt -o 16_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 24_start=6000_end=6250.fasta_exec.py
qsub -q TullerNano -r y  -e 17_job_error.txt -o 17_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 25_start=6250_end=6500.fasta_exec.py
qsub -q TullerNano -r y  -e 18_job_error.txt -o 18_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 26_start=6500_end=6750.fasta_exec.py
qsub -q TullerNano -r y  -e 19_job_error.txt -o 19_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 27_start=6750_end=7000.fasta_exec.py
qsub -q TullerNano -r y  -e 20_job_error.txt -o 20_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 28_start=7000_end=7250.fasta_exec.py
qsub -q TullerNano -r y  -e 21_job_error.txt -o 21_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 29_start=7250_end=7500.fasta_exec.py
qsub -q TullerNano -r y  -e 22_job_error.txt -o 22_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 2_start=500_end=750.fasta_exec.py
qsub -q TullerNano -r y  -e 23_job_error.txt -o 23_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 30_start=7500_end=7750.fasta_exec.py
qsub -q TullerNano -r y  -e 24_job_error.txt -o 24_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 31_start=7750_end=8000.fasta_exec.py
qsub -q TullerNano -r y  -e 25_job_error.txt -o 25_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 32_start=8000_end=8250.fasta_exec.py
qsub -q TullerNano -r y  -e 26_job_error.txt -o 26_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 33_start=8250_end=8500.fasta_exec.py
qsub -q TullerNano -r y  -e 27_job_error.txt -o 27_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 34_start=8500_end=8750.fasta_exec.py
qsub -q TullerNano -r y  -e 28_job_error.txt -o 28_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 35_start=8750_end=9000.fasta_exec.py
qsub -q TullerNano -r y  -e 29_job_error.txt -o 29_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 36_start=9000_end=9250.fasta_exec.py
qsub -q TullerNano -r y  -e 30_job_error.txt -o 30_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 37_start=9250_end=9500.fasta_exec.py
qsub -q TullerNano -r y  -e 31_job_error.txt -o 31_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 38_start=9500_end=9750.fasta_exec.py
qsub -q TullerNano -r y  -e 32_job_error.txt -o 32_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 39_start=9750_end=10000.fasta_exec.py
qsub -q TullerNano -r y  -e 33_job_error.txt -o 33_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 3_start=750_end=1000.fasta_exec.py
qsub -q TullerNano -r y  -e 34_job_error.txt -o 34_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 40_start=10000_end=10250.fasta_exec.py
qsub -q TullerNano -r y  -e 35_job_error.txt -o 35_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 41_start=10250_end=10500.fasta_exec.py
qsub -q TullerNano -r y  -e 36_job_error.txt -o 36_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 42_start=10500_end=10750.fasta_exec.py
qsub -q TullerNano -r y  -e 37_job_error.txt -o 37_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 43_start=10750_end=11000.fasta_exec.py
qsub -q TullerNano -r y  -e 38_job_error.txt -o 38_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 44_start=11000_end=11250.fasta_exec.py
qsub -q TullerNano -r y  -e 39_job_error.txt -o 39_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 45_start=11250_end=11500.fasta_exec.py
qsub -q TullerNano -r y  -e 40_job_error.txt -o 40_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 46_start=11500_end=11750.fasta_exec.py
qsub -q TullerNano -r y  -e 41_job_error.txt -o 41_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 47_start=11750_end=12000.fasta_exec.py
qsub -q TullerNano -r y  -e 42_job_error.txt -o 42_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 48_start=12000_end=12250.fasta_exec.py
qsub -q TullerNano -r y  -e 43_job_error.txt -o 43_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 49_start=12250_end=12456.fasta_exec.py
qsub -q TullerNano -r y  -e 44_job_error.txt -o 44_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 4_start=1000_end=1250.fasta_exec.py
qsub -q TullerNano -r y  -e 45_job_error.txt -o 45_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 5_start=1250_end=1500.fasta_exec.py
qsub -q TullerNano -r y  -e 46_job_error.txt -o 46_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 6_start=1500_end=1750.fasta_exec.py
qsub -q TullerNano -r y  -e 47_job_error.txt -o 47_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 7_start=1750_end=2000.fasta_exec.py
qsub -q TullerNano -r y  -e 48_job_error.txt -o 48_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 8_start=2000_end=2250.fasta_exec.py
qsub -q TullerNano -r y  -e 49_job_error.txt -o 49_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 9_start=2250_end=2500.fasta_exec.py
