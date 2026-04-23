#!/usr/bin/env python3
"""
HPC System Statistics Reporter

This utility script provides a quick snapshot of the current system resource
utilization on an HPC node. It reports CPU usage, memory availability, and
disk space information in a human-readable format.

This script is useful for:
    - Monitoring node availability before job submission
    - Diagnosing resource constraints
    - Checking disk space before running simulations
    - Validating node health in HPC environments

Dependencies:
    - psutil: Cross-platform library for retrieving system and process utilities
    
Usage:
    python3 get_system_stats.py
    
Output example:
    CPU Info-->  45.2%
    Memory Info--> 128.5MB free / 8192.0MB total ( 98.4% )
    Disk Info--> 250.3GB free / 1000.0GB total ( 75.0% )
"""

import psutil

# ==============================================================================
# CPU STATISTICS
# ==============================================================================

# Get current CPU utilization percentage across all cores
# Returns values from 0-100% representing average CPU usage since last call
cpu = str(psutil.cpu_percent()) + '%'

# ==============================================================================
# MEMORY STATISTICS
# ==============================================================================

# Retrieve virtual memory statistics for the system
# Returns a named tuple with: total, available, percent, used, free
memory = psutil.virtual_memory()

# Convert available memory from Bytes to MB
# Conversion: Bytes -> KB (÷1024) -> MB (÷1024)
available = round(memory.available / 1024.0 / 1024.0, 1)

# Convert total memory from Bytes to MB
total = round(memory.total / 1024.0 / 1024.0, 1)

# Format memory information string with available, total, and usage percentage
mem_info = str(available) + 'MB free / ' + str(total) + 'MB total ( ' + str(memory.percent) + '% )'

# ==============================================================================
# DISK STATISTICS
# ==============================================================================

# Retrieve disk usage statistics for the root filesystem
# The '/' parameter specifies the root directory (can be changed for other mounts)
disk = psutil.disk_usage('/')

# Convert free disk space from Bytes to GB
# Conversion: Bytes -> KB (÷1024) -> MB (÷1024) -> GB (÷1024)
free = round(disk.free / 1024.0 / 1024.0 / 1024.0, 1)

# Convert total disk space from Bytes to GB
total = round(disk.total / 1024.0 / 1024.0 / 1024.0, 1)

# Format disk information string with free space, total space, and usage percentage
disk_info = str(free) + 'GB free / ' + str(total) + 'GB total ( ' + str(disk.percent) + '% )'

# ==============================================================================
# OUTPUT SYSTEM STATISTICS
# ==============================================================================

print("CPU Info--> ", cpu)
print("Memory Info-->", mem_info)
print("Disk Info-->", disk_info)

