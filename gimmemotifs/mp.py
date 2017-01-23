from gimmemotifs.config import MotifConfig
import signal
from multiprocessing import Pool

n_cpus = int(MotifConfig().get_default_params()["ncpus"])

original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = Pool(processes=n_cpus, maxtasksperchild=1000)
signal.signal(signal.SIGINT, original_sigint_handler)
