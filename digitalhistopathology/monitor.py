#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import json
import os
import time
from threading import Thread

import GPUtil
import psutil
import torch


class Monitor(Thread):
    def __init__(self, delay=5):
        """Class to monitor CPU, RAM and GPU in a thread.

        Args:
            delay (int, optional): Every delay seconds the monitoring of the CPU, RAM and GPU is lauched. Defaults to 5.
        """
        super(Monitor, self).__init__()
        self.delay = delay
        self.running = True
        self.cpu_counts = psutil.cpu_count()
        self.gpu_available = torch.cuda.is_available()
        if self.gpu_available:
            print("GPU available")
        self.cpu_results = []
        self.gpu_results = []
        self.ram_results = []

    def cpu_usage(self):
        current_process = psutil.Process(os.getpid())
        cpu_perc = current_process.cpu_percent(interval=1) / self.cpu_counts
        self.cpu_results.append({"time": time.time(), "perc": cpu_perc})
        print("CPU utilization percentage: {}%".format(cpu_perc))

    def gpu_usage(self, index=0):
        gpus = GPUtil.getGPUs()
        print("Number of GPUs: {}".format(len(gpus)))
        gpu = gpus[index]
        gpu_perc = gpu.memoryUtil * 100
        gpu_used_gb = gpu.memoryUsed / 1000
        self.gpu_results.append({"time": time.time(), "perc": gpu_perc, "used_gb": gpu_used_gb})
        print("GPU utilization percentage: {}%".format(gpu_perc))

    def ram_usage(self):
        current_process = psutil.Process(os.getpid())
        memory_info = current_process.memory_full_info()
        # Convert bytes to gigabytes
        used_gb = memory_info.uss / (1024**3)
        total_gb = memory_info.vms / (1024**3)
        # Calculate the percentage of memory usage
        percentage_used = (memory_info.uss / memory_info.vms) * 100
        self.ram_results.append({"time": time.time(), "perc": percentage_used, "used_gb": used_gb})
        print("RAM utilization percentage: {}%".format(percentage_used))

    def run(self):
        while self.running:
            self.ram_usage()
            self.cpu_usage()
            if self.gpu_available:
                self.gpu_usage()
            time.sleep(self.delay)

    def stop(self):
        self.running = False

    def save_results(self, saving_folder):
        """Save the monitoring of CPU, RAM and GPU in 3 separated json file.

        Args:
            saving_folder (str): Path to the folder where to save the json result files.
        """
        with open(os.path.join(saving_folder, "cpu.json"), "w") as fp:
            json.dump(self.cpu_results, fp)
        with open(os.path.join(saving_folder, "ram.json"), "w") as fp:
            json.dump(self.ram_results, fp)
        if self.gpu_available:
            with open(os.path.join(saving_folder, "gpu.json"), "w") as fp:
                json.dump(self.gpu_results, fp)
