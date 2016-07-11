#!/usr/bin/env python

def basename(args):
    modified_args = args[:]
    for i in range(len(modified_args)):
        for tag in [ ("--",""),
                     ("static_Bz","sBz"),
                     ("c13_percentage","cp"),
                     ("max_cluster_size","mcs"),
                     ("scale_factor","sf"),
                     ("log10_samples","ls") ]:
            modified_args[i] = modified_args[i].replace(tag[0],tag[1])

    return "-".join(modified_args)
