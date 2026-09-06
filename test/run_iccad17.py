import os
import csv
import sys
import subprocess as sp

bench_list = [
  "des_perf_1",
  "des_perf_a_md1",
  "des_perf_a_md2",
  "des_perf_b_md1",
  "des_perf_b_md2",
  "edit_dist_1_md1",
  "edit_dist_a_md2",
  "edit_dist_a_md3",
  "fft_2_md2",
  "fft_a_md2",
  "fft_a_md3",
  "pci_bridge32_a_md1",
  "pci_bridge32_a_md2",
  "pci_bridge32_b_md1",
  "pci_bridge32_b_md2",
  "pci_bridge32_b_md3"
]

def get_data(file):

  S = 0
  S_am = 0
  S_hpwl = 0
  M_max = 0
  overlap = 0
  N_p = 0
  N_e = 0
  runtime = 0
  avg_disp_sites = 0

  with open(file, 'r') as file_open:
    file_lines = file_open.readlines()

    for line in file_lines:
      line = line.strip()
      items = line.split()

      if len(items) == 0:
        continue

      if items[0] == "S_am":
        S_am = float(items[2])

      if items[0] == "Overlap":
        overlap = int(items[2])

      if items[0] == "Delta":
        S_hpwl = float(items[3])

      if items[0] == "N_p":
        N_p = int(items[2])

      if items[0] == "N_e":
        N_e = int(items[2])

      if items[0] == "Max" and items[1] == "Displacement" and items[4] == "(rows)":
        M_max = float(items[3])

      if items[0] == "Avg" and items[1] == "Displacement" and items[4] == "(sites)":
        avg_disp_sites = float(items[3])

      if items[0] == "Legalization" and items[1] == "Time":
        runtime = float(items[3])

      if items[0] == "S":
        S = float(items[2])

  return [S, S_am, S_hpwl, M_max, N_p, N_e, avg_disp_sites, overlap, runtime]

def write_tcl(design_name, template_name, tcl_name):
  inFileName = template_name
  inFile  = open(inFileName, "r")

  outFileName = tcl_name
  outFile = open(outFileName, "w")

  text = inFile.read()
  text = text.replace("__DESIGN_NAME__", design_name)

  outFile.write(text)
  outFile.close()
  inFile.close()

def run_and_eval(design_name):

  template_name = "template_iccad17.tcl"
  tcl_name = "temp_tcl.tcl"

  write_tcl(design_name, template_name, tcl_name)

  log_name = "temp_log.txt"
  run_cmd  = f"./SkyLine {tcl_name} | tee {log_name}"

  sp.call(run_cmd, shell=True, stdout=sp.DEVNULL)

  S, S_am, S_hpwl, M_max, N_p, N_e, avg_disp_sites, overlap, runtime = get_data(log_name)

  print(f"{design_name}")
  print(f"S: {S}")
  print(f"S_am: {S_am}")
  print(f"S_hpwl: {S_hpwl}")
  print(f"M_max: {M_max}")
  print(f"N_p: {N_p}")
  print(f"N_e: {N_e}")
  print(f"Avg Disp: {avg_disp_sites} (sites)")
  print(f"Overlap: {overlap}")
  print(f"Runtime: {runtime}")

  sp.call(f"rm {log_name}", shell=True)
  sp.call(f"rm {tcl_name}", shell=True)

  return [design_name, S, S_am, S_hpwl, M_max, N_p, N_e, avg_disp_sites, overlap, runtime]

csv_name = "iccad17_results.csv"
if os.path.exists(csv_name):
  sp.call(f"rm {csv_name}", shell=True)

with open(csv_name, "a") as csv_file:
  csv_wr = csv.writer(csv_file)
  csv_wr.writerow(["Benchmark", "S", "S_am", "S_hpwl", "M_max", "N_p", "N_e", "AvgDisp (sites)", "Overlap", "Runtime"])

for design in bench_list:
  data = run_and_eval(design)
  with open(csv_name, "a") as csv_file:
    csv_wr = csv.writer(csv_file)
    csv_wr.writerow(data)
