import os
from utilis import seed_it
configs = ['config_rf', 'config_svr', 'config_dt', 'config_nn', 'config_lasso', 'config_ridge', 'config_lr']
seeds = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

for config in configs:
    for seed in seeds:
        seed_it(seed)
        os.system('python main_sup.py --config configs/{}.json'.format(config))
        break
    break