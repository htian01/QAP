import os
import subprocess

load_path = "./dataset/qaplib/problems_dat"
save_path = "./dataset/qaplib/hg98_reparam_it2000"
log_path = "./dataset/qaplib/hg98_reparam_it2000"

def run_problem(problem_name):
    print(f"Run {problem_name}")
    cmd = f"nohup ./main {load_path} {save_path} {problem_name} 2000"
    with open(os.path.join(log_path, problem_name + "_reparam.log"), 'w') as f:
        subprocess.call(cmd.split(), stdout=f)

logfiles = set(file[:-len("_reparam.dd")] for file in os.listdir(save_path) if file.endswith(".dd"))

for filename in os.listdir(load_path):
    problem_name = filename[:-4]
    n = int("".join(list(filter(str.isdigit, problem_name))))
    if n < 200:
        if problem_name not in logfiles:
            run_problem(problem_name)



