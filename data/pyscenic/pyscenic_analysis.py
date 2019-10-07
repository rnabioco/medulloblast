

import os
import sys 
import glob
import pickle
import pandas as pd
import numpy as np


from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

if __name__ == '__main__':
  #cluster = LocalCluster(n_workers = 44, threads_per_worker = 1,
  #        memory_limit=8e9)

  #client = Client(cluster)
  
  #SCHEDULER=client.scheduler.address
  DATABASE_DIR = "../../dbases/pyscenic" 
  TFS_FNAME = os.path.join(DATABASE_DIR, 'tfs.txt')
  RESOURCES_DIR = "tables"

  # download tf list
  import urllib.request
  
  if not os.path.isfile(TFS_FNAME): 
    tf_url = "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_curated_tfs.txt"
    urllib.request.urlretrieve(tf_url, TFS_FNAME)

  DATA_DIR="output_grp34"
  DATABASES_GLOB = os.path.join(DATABASE_DIR, "hg38_*.mc9nr.feather")
  MOTIF_ANNOTATIONS_FNAME = os.path.join(DATABASE_DIR,
          "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
  SC_EXP_FNAME = os.path.join(RESOURCES_DIR, "gp34_expr_matrix.tsv.gz") 
  REGULONS_FNAME = os.path.join(DATA_DIR, "regulons.p")
  MOTIFS_FNAME = os.path.join(DATA_DIR, "motifs.csv")
  AUC_RES = os.path.join(DATA_DIR, "auc_res.csv")
  
  print("loading matrix", file=sys.stdout)
  ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
  ex_matrix.shape

  #tf_names = load_tf_names(TFS_FNAME)

  #db_fnames = glob.glob(DATABASES_GLOB)
  #def name(fname):
  #    return os.path.splitext(os.path.basename(fname))[0]
  #dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
  #dbs

  print("loaded dbs", file=sys.stdout)
  #print(dbs, file=sys.stdout)
  print("running grnboost", file=sys.stdout)
  #adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True,
  #        client_or_address=client)
  #
  #adj_tmp_file = os.path.join(DATA_DIR, "adj.pkl")

  #with open(adj_tmp_file, "wb") as f:
  #    pickle.dump(adjacencies, f)

  #modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
  #
  #with open(os.path.join(DATA_DIR, "modules.pkl"), "wb") as f:
  #    pickle.dump(modules, f)
  
  #with open(os.path.join(DATA_DIR, "modules.pkl"), "rb") as f:
  #    modules = pickle.load(f)
  
  print("running prune2df", file=sys.stdout)
  #df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, nes_threshold=2.5,
  #        client_or_address=client)
  
  #df.to_csv(os.path.join(DATA_DIR, "tmp_df.csv"))
  print("running df2regulons", file=sys.stdout)
  # Create regulons from this table of enriched motifs.
  #regulons = df2regulons(df)
  
  # Save the enriched motifs and the discovered regulons to disk.
  #df.to_csv(MOTIFS_FNAME)
  #with open(REGULONS_FNAME, "wb") as f:
  #    pickle.dump(regulons, f)
  
  with open(REGULONS_FNAME, "rb") as f:
      regulons = pickle.load(f)
  
  #client.close()
  #cluster.close()
  print("running aucell", file=sys.stdout)
  auc_mtx = aucell(ex_matrix, regulons, num_workers=11)
  
  auc_mtx.to_csv(AUC_RES)
