import pandas as pd
import matplotlib.pyplot as plt

import yaml
import os

with open('conf.yaml') as f:
    
    configs = yaml.load(f, Loader=yaml.FullLoader)

    for routing_type in configs['custom_routing_type']:
        if routing_type == "NONE":
            for routing in configs['routing']:
                if routing == "DYAD":
                    cmd = "~/noxim/noxim_with_tasks/bin/noxim -config ../configs/{}.yaml -power ~/noxim/noxim_with_tasks/bin/power.yaml -routing {} 0.6 -seed 7688787787 \
                         > out_file.txt".format("mcsl_default", routing)
                    os.system(cmd)
                    continue
                cmd = "~/noxim/noxim_with_tasks/bin/noxim -config ../configs/{}.yaml -power ~/noxim/noxim_with_tasks/bin/power.yaml -routing {} -seed 7688787787 \
                     > out_file.txt".format("mcsl_default", routing)
                print(cmd)
                os.system(cmd)
        else:
            cmd = "~/noxim/noxim_with_tasks/bin/noxim -config ../configs/{}.yaml -power ~/noxim/noxim_with_tasks/bin/power.yaml -custom {} -seed 7688787787 \
                 > out_file.txt".format("mcsl_default", routing_type)
            os.system(cmd)

