#!/usr/bin/env python3
import sys, zipfile, os, shutil, subprocess
from main import workflow
in_file = sys.argv[1]          # 用户上传的文件
out_dir = sys.argv[2]          # 结果输出目录

# 伪处理：生成 3 个结果文件
os.makedirs(out_dir, exist_ok=True)
# for i in 1,2,3:
#     with open(f'{out_dir}/result_{i}.txt','w') as f:
#         f.write(f'这是第 {i} 个结果，原文件={in_file}\n')
workflow(name='fastqc', manual=in_file, run_help_command=False, output_location=out_dir)
