import threading
import time
from contextlib import contextmanager

import pandas as pd
import psutil


class Benchmark:
    """
    Records wall time, CPU time and peak memory of named steps.

    Memory and CPU time are summed over this process and all its child processes (e.g. the
    workers of a dask LocalCluster). They are sampled every `interval` seconds, so peaks
    shorter than that and children that start and end between two samples are missed.

    Usage:
        bench = Benchmark()
        with bench.step('inference', n_targets=100):
            ...
        bench.table().to_csv('benchmark.tsv', sep='\t', index=False)
    """

    def __init__(self, interval=0.5):
        self.interval = interval
        self.process = psutil.Process()
        self.records = []

    def _processes(self):
        return [self.process] + self.process.children(recursive=True)

    def _sample(self, cpu, peak):
        rss = 0
        for p in self._processes():
            try:
                with p.oneshot():
                    rss += p.memory_info().rss
                    times = p.cpu_times()
                    cpu[p.pid] = times.user + times.system
            except psutil.NoSuchProcess:
                pass
        peak[0] = max(peak[0], rss)

    @contextmanager
    def step(self, name, **info):
        """
        Measures the code in the with block as step `name`. Keyword arguments are added as
        columns (e.g. input sizes). More can be added inside the block via the yielded dict.
        """
        start_cpu, cpu, peak = {}, {}, [0]
        self._sample(start_cpu, peak)
        stop = threading.Event()

        def sample():
            while not stop.wait(self.interval):
                self._sample(cpu, peak)

        thread = threading.Thread(target=sample, daemon=True)
        start = time.perf_counter()
        thread.start()
        try:
            yield info
        finally:
            stop.set()
            thread.join()
            self._sample(cpu, peak)
            wall = time.perf_counter() - start
            cpu_time = sum(t - start_cpu.get(pid, 0) for pid, t in cpu.items())
            self.records.append({'step': name, 'wall_time_s': wall, 'cpu_time_s': cpu_time,
                                 'peak_rss_gb': peak[0] / 1e9, 'n_processes': len(cpu), **info})
            print(f"[benchmark] {name}: {wall:.1f} s wall, {cpu_time:.1f} s CPU, {peak[0] / 1e9:.2f} GB peak RSS")

    def table(self):
        """
        Returns:
            pd.DataFrame: One row per step with 'step', 'wall_time_s', 'cpu_time_s',
                'peak_rss_gb', 'n_processes' and the info columns, plus a 'total' row.
        """
        steps = pd.DataFrame(self.records)
        total = {'step': 'total', 'wall_time_s': steps['wall_time_s'].sum(),
                 'cpu_time_s': steps['cpu_time_s'].sum(), 'peak_rss_gb': steps['peak_rss_gb'].max(),
                 'n_processes': steps['n_processes'].max()}
        return pd.concat([steps, pd.DataFrame([total])], ignore_index=True)
