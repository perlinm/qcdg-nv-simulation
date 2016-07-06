#!/usr/bin/env python

def basename(args):
    base_args = args[:6]
    opts = args[6:]
    return ("{}-sBz-{}-cp-{}-mcs-{}-sf-{}-ls-{}".format(*base_args)
            + "".join([ "-"+o for o in opts ]))
