import os
from utilis import seed_it
configs = ['mean', 'med',  'geo', 'l2', 'mc1', 'mc2', 'mc3','mc4','thurstone','stuart', 'rra', 'bard', 'birra', 'big', 'cemc_kendall', 'cemc_spearman']
seeds = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
path = '/mnt/c/Users/gqu/OneDrive - UTHealth Houston/projects/Genevic/data/data_dict.pkl'
for config in configs:
    # for seed in seeds:
    #     seed_it(seed)
    os.system('python main_unsup.py --data_path "{}" --config {} --prefix {}'.format(path, config, config))