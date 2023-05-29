from ogtk.utils import db
import matplotlib.pyplot as plt
import seaborn as sns
import ogtk.utils as ut
import ogtk.ltr.shltr as shltr
import numpy as np
from typing import Sequence,Optional,Iterable
import polars as pl
import anndata as ad
import os
from functools import wraps


class Xp(db.Xp):
    def hola(self):
        print("hola")

def rotate_labs(fg, ylim=None):
    if len(fg.axes_dict)>0:
        for ax in fg.axes_dict.values():
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.grid()
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set_axisbelow(True)
    else:
        ax = fg.ax
        ax.grid()
        ax.set_axisbelow(True)
        art = fg.ax.set_xticklabels(fg.ax.get_xticklabels(), rotation=90)
        if ylim is not None:
            ax.set_ylim(ylim)
