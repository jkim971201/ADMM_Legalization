import os
import csv
import sys
import subprocess as sp

bench_list = [
  "mgc_superblue11_a",
  "mgc_superblue12",
  "mgc_superblue14",
  "mgc_superblue16_a",
  "mgc_superblue19"
]

def get_data(file):

  total_runtime = 0
  preplace_runtime = 0
  admm_runtime = 0
  post_runtime = 0
  max_disp = 0
  total_disp = 0

  with open(file, 'r') as file_open:
    file_lines = file_open.readlines()

    for line in file_lines:
      line = line.strip()
      items = line.split()

      if len(items) == 0:
        continue

      if items[0] == "Sum" and items[4] == "(sites)":
        total_disp = float(items[3])
      elif items[0] == "Max" and items[4] == "um":
        max_disp = float(items[3])
      elif items[0] == "preplace":
        preplace_runtime = float(items[3])
      elif items[0] == "admmLegalize":
        admm_runtime = float(items[3])
      elif items[0] == "postProcess":
        post_runtime = float(items[3])

  total_runtime = preplace_runtime + admm_runtime + post_runtime
  return [total_disp, max_disp, total_runtime]

def write_tcl(design_name, template_name, tcl_name):
  inFileName = template_name
  inFile  = open(inFileName, "r")

  outFileName = tcl_name
  outFile = open(outFileName, "w")

  text = inFile.read()
  text = text.replace("__DESIGN__", design_name)

  outFile.write(text)
  outFile.close()
  inFile.close()

def run_and_eval(design_name):

  template_name = "__skyline_ispd15_template.tcl"
  tcl_name = "temp_tcl.tcl"

  write_tcl(design_name, template_name, tcl_name)

  log_name = "temp_log.txt"
  run_cmd  = f"./SkyLine {tcl_name} | tee {log_name}"

  #sp.call(run_cmd, shell=True, stdout=sp.DEVNULL)
  sp.call(run_cmd, shell=True)

  total_disp, max_disp, runtime = get_data(log_name)

  print(f"{design_name}")
  print(f"TotalDisp : {total_disp} sites")
  print(f"Max  Disp : {max_disp} um")
  print(f"Runtime   : {runtime} s")

  sp.call(f"rm {log_name}", shell=True)
  sp.call(f"rm {tcl_name}", shell=True)

  return [design_name, total_disp, max_disp, runtime]

csv_name = "ispd15_results.csv"
if os.path.exists(csv_name):
  sp.call(f"rm {csv_name}", shell=True)

with open(csv_name, "a") as csv_file:
  csv_wr = csv.writer(csv_file)
  csv_wr.writerow(["Benchmark", "TotalDisp (sites)", "Maxdisp (um)", "Runtime"])

for design in bench_list:
  data = run_and_eval(design)
  with open(csv_name, "a") as csv_file:
    csv_wr = csv.writer(csv_file)
    csv_wr.writerow(data)
