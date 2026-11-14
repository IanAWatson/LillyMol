#!/usr/bin/env bash

# Write to stdout a suggested value for THREADS: how many parallel
# compile jobs (make -j, bazel --jobs) the current machine can sustain.
# Diagnostic messages must go to stderr because callers capture stdout.

# Logical core count. nproc respects cgroup cpusets, so this is already
# correct under most container CPU-pinning setups; it does not account
# for a cgroup v2 cpu.max *quota* (as opposed to a cpuset), so check
# that too when present.
cores=$(nproc 2>/dev/null || echo 1)

if [[ -r /sys/fs/cgroup/cpu.max ]] ; then
    # cgroup v2: "<quota> <period>" in microseconds, or "max <period>".
    read -r cpu_quota cpu_period < /sys/fs/cgroup/cpu.max
    if [[ "${cpu_quota}" != 'max' && "${cpu_period}" -gt 0 ]] ; then
        quota_cores=$(( cpu_quota / cpu_period ))
        if [[ ${quota_cores} -ge 1 && ${quota_cores} -lt ${cores} ]] ; then
            cores=${quota_cores}
        fi
    fi
elif [[ -r /sys/fs/cgroup/cpu/cpu.cfs_quota_us && -r /sys/fs/cgroup/cpu/cpu.cfs_period_us ]] ; then
    # cgroup v1 equivalent.
    cpu_quota=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us)
    cpu_period=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)
    if [[ ${cpu_quota} -gt 0 && ${cpu_period} -gt 0 ]] ; then
        quota_cores=$(( cpu_quota / cpu_period ))
        if [[ ${quota_cores} -ge 1 && ${quota_cores} -lt ${cores} ]] ; then
            cores=${quota_cores}
        fi
    fi
fi

# Heavily templated code (eg anything that pulls in Eigen) can need well
# over 1GB per compile job. Cap parallelism by available memory too, so
# a box with many cores but little RAM does not thrash/OOM instead of
# just going a bit slower. Budget ~2GB/job.
mem_kb=$(awk '/MemAvailable:/ {print $2}' /proc/meminfo 2>/dev/null)
if [[ -z "${mem_kb}" ]] ; then
    mem_kb=$(( 2 * 1024 * 1024 )) # unknown: assume 2GB, forces jobs_by_mem=1
fi
if [[ -r /sys/fs/cgroup/memory.max ]] ; then
    # cgroup v2
    cgroup_mem=$(cat /sys/fs/cgroup/memory.max)
    if [[ "${cgroup_mem}" != 'max' ]] ; then
        cgroup_mem_kb=$(( cgroup_mem / 1024 ))
        if [[ ${cgroup_mem_kb} -lt ${mem_kb} ]] ; then
            mem_kb=${cgroup_mem_kb}
        fi
    fi
elif [[ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]] ; then
    # cgroup v1. A very large value here means "no limit set".
    cgroup_mem=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes)
    if [[ ${cgroup_mem} -lt $(( 1024 * 1024 * 1024 * 1024 )) ]] ; then
        cgroup_mem_kb=$(( cgroup_mem / 1024 ))
        if [[ ${cgroup_mem_kb} -lt ${mem_kb} ]] ; then
            mem_kb=${cgroup_mem_kb}
        fi
    fi
fi

jobs_by_mem=$(( mem_kb / (2 * 1024 * 1024) ))
if [[ ${jobs_by_mem} -lt 1 ]] ; then
    jobs_by_mem=1
fi

threads=${cores}
if [[ ${jobs_by_mem} -lt ${threads} ]] ; then
    threads=${jobs_by_mem}
fi

echo "cores=${cores} jobs_by_mem=${jobs_by_mem} -> THREADS=${threads}" >&2
echo "${threads}"
