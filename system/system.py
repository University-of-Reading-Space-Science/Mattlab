#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 17:09:17 2024

System functions

@author: vy902033
"""

import os as os



def _setup_dirs_():
    
    cwd = os.path.abspath(os.path.dirname(__file__))
    root = os.path.dirname(cwd)
    datapath = os.path.join(root,'Data')
    
    dirs = {"datapath": datapath}
    dirs['root'] = root
    
    return dirs