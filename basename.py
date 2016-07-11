#!/usr/bin/env python

log10_samples_hook = "log10_samples"

def basename(args):
    modified_args = [ str(a) for a in args[:] ]
    for i in range(len(modified_args)):
        for tag in [ ("--",""),
                     ("static_Bz","sBz"),
                     ("c13_percentage","cp"),
                     ("max_cluster_size","mcs"),
                     ("scale_factor","sf"),
                     (log10_samples_hook,"ls") ]:
            modified_args[i] = modified_args[i].replace(tag[0],tag[1])

    return "-".join(modified_args)
