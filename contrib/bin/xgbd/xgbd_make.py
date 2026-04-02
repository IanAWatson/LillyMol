# Build and commit an xgboost model
# Deliberately simplistic in approach

import os
import subprocess
import sys
from typing import List

import numpy as np
import optuna
from optuna.samplers import TPESampler
import pandas as pd

from absl import app, flags, logging
from google.protobuf import text_format
from google.protobuf import json_format
from matplotlib import pyplot
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from xgboost import XGBClassifier, XGBRegressor, plot_importance
# import xgbd.xgboost_model_pb2
from xgbd import xgboost_model_pb2
from xgbd import class_label_translation_pb2

FLAGS = flags.FLAGS

flags.DEFINE_integer("min_points", 10, "do NOT build a model if there are fewer than min_points")
flags.DEFINE_string("activity", "", "Name of training set activity file")
flags.DEFINE_boolean("classification", False, "True if this is a classification task")
flags.DEFINE_string("mdir", "", "Directory into which the model is placed")
flags.DEFINE_integer("max_num_features", 0, "Maximum number of features to plot in variable importance")
flags.DEFINE_string("feature_importance", "", "Compute feature importance. Use 'def' to use default file")
flags.DEFINE_boolean("optuna", False, "File containing optuna config")
flags.DEFINE_integer("xgverbosity", 0, "xgboost verbosity")
flags.DEFINE_string("proto", "", "A file containing an XGBoostParameters proto")
flags.DEFINE_float("eta", 0.4, "xgboost learning rate parameter eta")
flags.DEFINE_integer("max_depth", 5, "xgboost max depth")
flags.DEFINE_integer("n_estimators", 500, "xboost number of estimators")
flags.DEFINE_float("subsample", 1.0, "subsample ratio for training instances")
flags.DEFINE_float("min_child_weight", 1.0, "subsample ratio for training instances")
flags.DEFINE_float("colsample_bytree", 1.0, "subsampling occurs once for every tree constructed")
flags.DEFINE_float("colsample_bylevel", 1.0, "subsampling occurs once for every new depth level reached")
flags.DEFINE_float("colsample_bynode", 1.0, "subsampling occurs once for every time a new split is evaluated")
flags.DEFINE_float("reg_lambda", 1.0, "L1 regularization")
flags.DEFINE_float("reg_alpha", 0.0, "L2 regularization")
flags.DEFINE_float("gamma", 0.0, "minimum loss reduction")
flags.DEFINE_enum("tree_method", "auto", ["auto", "exact", "approx", "hist"], "tree construction method: auto exact approx hist")
flags.DEFINE_integer("nthreads", 8, "number of threads to use, default is 8")
flags.DEFINE_boolean("rescore", False, "Rescore the training set to establish linear correction function")


class Options:
  def __init__(self):
    self.min_points = 10
    self.classification = False
    self.mdir: str = ""
    self.max_num_features: int = 15
    self.descriptor_fname: str = ""
    self.activity_fname: str = ""
    self.verbosity = 0
    self.proto = xgboost_model_pb2.XGBoostParameters()
    self.optuna = ""

  def read_proto(self, fname)->bool:
    """Read self.proto from `fname`
    """
    with open(fname, "r") as reader:
      text = reader.read()

    self.proto = text_format.Parse(text, xgboost_model_pb2.XGBoostParameters())
    if not self.proto:
      logging.error("Cannot intpret %s", text)
      return False

    return True

def to_array(input:str) -> List[str]:
  return input.split(' ')

def write_class_label_translation(options: Options, categories:np.array, class_counts:np.array)->bool:
  """Write the class label proto to the model directory.
  """
  proto = class_label_translation_pb2.ClassLabelTranslation()
  proto.to_numeric[categories[0]] = 0
  proto.to_numeric[categories[1]] = 1


  fname = os.path.join(options.mdir, "class_label_translation.dat")
  with open(fname, "wb") as output:
    serialised = proto.SerializeToString()
    output.write(serialised)

  fname = os.path.join(options.mdir, "class_label_translation.json")
  with open(fname, "w") as output:
    output.write(json_format.MessageToJson(proto))

  return True

def classification(x, y, options: Options)->bool:
  """build a classification model
    Args:
      x: feature matrix
      y: response - must be translated to 0,1. Not implemented...
  """
  categories, counts = np.unique(y, return_counts=True)
  print(type(y))
  if len(categories) != 2:
    logging.error("Must be two classes %d not possible\n", len(categories))
    return False
  if counts[0] <= counts[1]:
    logging.info("less")
    y = (y == categories[0]).astype(int)
    print(y)
  else:
    y = (y != categories[0]).astype(int)
    counts[0], counts[1] = counts[1], counts[0]
    categories[0], categories[1] = categories[1], categories[0]
    logging.info("greater")
  write_class_label_translation(options, categories, counts)
  print(categories)
  print(counts)
  booster = XGBClassifier(verbosity=options.verbosity)
  booster.fit(x, y)

  booster.save_model(os.path.join(options.mdir, "xgboost.json"))

  return True

# Optuna objective function
def objective(trial):
    param = {
        'booster': trial.suggest_categorical('booster', ['gbtree', 'dart']),
        'lambda': trial.suggest_float('lambda', 1e-8, 1.0, log=False),
        'max_depth': trial.suggest_int('max_depth', 3, 9),
        'eta': trial.suggest_float('eta', 1e-8, 1.0, log=True),
        nthreads: nthreads
    }
    model = xgb.XGBClassifier(**param)
    return cross_val_score(model, X_train, y_train, cv=3).mean()

def regression_with_optuna(x, y, options: Options):
  """build a regression model.
    Thanks ChatGPT
  """
  match options.proto.tree_method:
    case xgboost_model_pb2.AUTO:
      tree_method = 'auto'
    case xgboost_model_pb2.EXACT:
      tree_method = 'exact'
    case xgboost_model_pb2.APPROX:
      tree_method = 'approx'
    case xgboost_model_pb2.HIST:
      tree_method = 'hist'
    case _:
      tree_method = 'hist'

  x = np.asarray(x)
  y = np.asarray(y)
  nthread = None
  if FLAGS.nthreads > 0:
    nthread = FLAGS.nthreads

  n_splits = 5

  random_state = None
  kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

  def objective(trial: optuna.Trial) -> float:
      params = {
          # Core
          "n_estimators": 10_000,  # large; early stopping finds the effective number
          "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
          "max_depth": trial.suggest_int("max_depth", 2, 12),
          "min_child_weight": trial.suggest_float("min_child_weight", 1e-3, 50.0, log=True),
          "subsample": trial.suggest_float("subsample", 0.5, 1.0),
          "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
          "gamma": trial.suggest_float("gamma", 0.0, 10.0),
          "reg_alpha": trial.suggest_float("reg_alpha", 1e-10, 10.0, log=True),
          "reg_lambda": trial.suggest_float("reg_lambda", 1e-10, 100.0, log=True),

          # Tree method (change if you want GPU)
          "tree_method": "hist",

          # Objective / metric
          "objective": "reg:squarederror",
          "eval_metric": "rmse",

          # Repro / speed
          "random_state": random_state,
          "n_jobs": nthread,
          "early_stopping_rounds": 200,
      }

      fold_rmses = []
      for fold_idx, (tr_idx, va_idx) in enumerate(kf.split(x), start=1):
          x_tr, x_va = x[tr_idx], x[va_idx]
          y_tr, y_va = y[tr_idx], y[va_idx]

          model = XGBRegressor(**params)

          # Early stopping: use a validation set within each CV fold
          model.fit(
              x_tr,
              y_tr,
              eval_set=[(x_va, y_va)],
              verbose=False,
          )

          preds = model.predict(x_va)
          rmse = mean_squared_error(y_va, preds)
          fold_rmses.append(rmse)

          # Let Optuna prune unpromising trials
          trial.report(float(np.mean(fold_rmses)), step=fold_idx)
          if trial.should_prune():
              raise optuna.TrialPruned()

      return float(np.mean(fold_rmses))

  study = optuna.create_study(
      direction="minimize",
      sampler=TPESampler(seed=random_state),
      pruner=optuna.pruners.MedianPruner(n_startup_trials=10, n_warmup_steps=1),
  )

  n_trials = 180
  timeout = 6000
  study.optimize(objective, n_trials=n_trials, timeout=timeout, show_progress_bar=False)

  print(study.best_params)
  return study

 
def regression(x, y, options: Options):
  """build a regression model.
  """
  match options.proto.tree_method:
    case xgboost_model_pb2.AUTO:
      tree_method = 'auto'
    case xgboost_model_pb2.EXACT:
      tree_method = 'exact'
    case xgboost_model_pb2.APPROX:
      tree_method = 'approx'
    case xgboost_model_pb2.HIST:
      tree_method = 'hist'
    case _:
      tree_method = 'hist'

  if options.optuna:
    return regression_with_optuna(x, y, options)
  
  nthread = None
  if FLAGS.nthreads > 0:
    nthread = FLAGS.nthreads

  booster = XGBRegressor(verbosity=options.verbosity,
                eta=options.proto.eta,
                max_depth=options.proto.max_depth,
                n_estimators=options.proto.n_estimators,
                colsample_bytree=options.proto.colsample_bytree,
                colsample_bylevel=options.proto.colsample_bylevel,
                colsample_bynode=options.proto.colsample_bynode,
                subsample=options.proto.subsample,
                reg_alpha=options.proto.reg_alpha,
                reg_lambda=options.proto.reg_lambda,
                gamma=options.proto.gamma,
                tree_method=tree_method,
                nthread=nthread
                )

  booster.fit(x, y)

  booster.save_model(os.path.join(options.mdir, "xgboost.json"))
  logging.info("Saved model to %s", os.path.join(options.mdir, "xgboost.json"))

  if options.max_num_features:
    plot_importance(booster, max_num_features=options.max_num_features)
    pyplot.show()

  if len(options.feature_importance) > 0:
    for itype in ["weight", "gain", "cover"]:
      feature_importance = booster.get_booster().get_score(importance_type=itype)
      feature_importance = sorted(feature_importance.items(), key=lambda x:x[1], reverse=True)

      fname = options.feature_importance
      if fname == "def" or fname == "DEF":
        fname = os.path.join(options.mdir, f'feature_importance.{itype}.txt')
      else:
        # Or should this be placed in `mdir` by default?
        fname = f"{fname}.{itype}.txt"

      with open(fname, "w") as writer:
        print("Feature Weight", file=writer)
        for f, i in feature_importance:
          print(f"{f} {i}", file=writer)

  return True


def build_xgboost_model(descriptor_fname: str,
                        activity_fname: str,
                        options: Options)->bool:
  """Build an xgboost model on the data in `descriptor_fname` and
     `activity_fname`.
    This function does data preprocessing.
  """

  descriptors = pd.read_csv(descriptor_fname, sep=' ', header=0, low_memory=False, na_values=['.'])
  logging.info("Read %d rows and %d columns from %s", len(descriptors),
                descriptors.shape[1], descriptor_fname)
  activity = pd.read_csv(activity_fname, sep=' ', header=0)
  logging.info("Read %d rows from %s", activity.shape[0], activity_fname)


  descriptors.rename(columns={descriptors.columns[0]: "Name"}, inplace=True)
  activity.rename(columns={activity.columns[0]: "Name"}, inplace=True)
# combined = pd.concat([activity.set_index("Name"),
#                       descriptors.set_index("Name")], axis=1, join='inner').reset_index() 
  combined = pd.merge(activity.set_index("Name"),
                        descriptors.set_index("Name"), how='inner', on=["Name"]).reset_index() 
  if len(combined) != len(descriptors):
    logging.error("Combined set has %d rows, need %d", len(combined), len(descriptors))
#   return False

  if len(combined) < options.min_points:
    logging.info("Not enough rows in training set %d", len(combined))
    return True  # Or should this be False? This is an OK exit, just no model is built

  if not os.path.isdir(options.mdir):
    os.mkdir(options.mdir)

  combined.to_csv(os.path.join(options.mdir, "train.xy"), sep= ' ', index=False)

  y = combined.iloc[:,1].to_numpy()

  x = combined.iloc[:,2:]
  features = x.columns
  x.apply(pd.to_numeric).to_numpy()

  options.descriptor_fname = descriptor_fname
  options.activity_fname = activity_fname

  rc = False
  if options.classification:
    rc = classification(x, y, options)
  else:
    rc = regression(x, y, options)

  if not rc:
    logging.info("Model did not build")
    return False

  response = activity.columns[1]

  proto = xgboost_model_pb2.XGBoostModel()
  proto.model_type = "XGBD"
  proto.classification = options.classification
  proto.response = response
  proto.parameters.CopyFrom(options.proto)

  for (column, feature) in enumerate(features):
    proto.name_to_col[feature] = column

  with open(os.path.join(options.mdir, "model_metadata.txt"), "w") as f:
    f.write(text_format.MessageToString(proto))
  with open(os.path.join(options.mdir, "model_metadata.dat"), "wb") as f:
    f.write(proto.SerializeToString())

  return True

def option_present(flag)->bool:
  """Return true if the option `flag` is in sys.argv
    Args:
      argv: usually the command line
      flag: a command line option. We look for -flag and --flag in argv.
  """
  if '-'+flag in sys.argv:
    return True
  if '--'+flag in sys.argv:
    return True
  return False


def rescore_training_set(options)->bool:
  """ A model has just been build in `options.mdir`.
      Rescore the training set and store the results.
  """
  train_pred = os.path.join(options.mdir, 'train.pred')
  with open(train_pred, 'w') as output:
    cmd = f"xgbd_evaluate.sh -mdir {options.mdir} {options.descriptor_fname}"
    subprocess.run(to_array(cmd), stdout=output, text=True)

  if not os.path.exists(train_pred):
    logging.error("%s did not create %s", cmd, train_pred)
    return False

  train_stats = os.path.join(options.mdir, 'train.stats')
  rescaling = os.path.join(options.mdir, 'rescaling.textproto')
  cmd = f"iwstats.sh -Y allequals -w -E {options.activity_fname} -p 2 -C {rescaling} {train_pred}"
  with open(train_stats, 'w') as output:
    subprocess.run(to_array(cmd), stdout=output, text=True)

  if not os.path.exists(rescaling):
    logging.error("%s did not create %s", cmd, rescaling)
    return False

  return True

def main(argv):
  """Build xgboost models from activity file and descriptor file.
  """
  if not FLAGS.activity:
    logging.error("Must specifythe name of the activity file with the --activity option")
    return False
  if len(argv) == 1:
    logging.error("Must specifythe name of the descriptor file as argument")
    return False
  if not FLAGS.mdir:
    logging.error("Must specify the model directory via the --mdir option")
    return False

  options = Options()
  options.min_points = FLAGS.min_points
  options.classification = FLAGS.classification
  options.mdir = FLAGS.mdir
  options.max_num_features = FLAGS.max_num_features
  options.feature_importance = FLAGS.feature_importance
  options.optuna = FLAGS.optuna
  options.verbosity = FLAGS.xgverbosity

  if FLAGS.classification and FLAGS.rescore:
    logging.error("The --classification and --rescore options are incompatible")
    return False

  # Build the proto first.
  # After that is done, we check for command line arguments that would
  # over-ride what has come in from the proto.
  if FLAGS.proto:
    if not options.read_proto(FLAGS.proto):
      logging.error("Cannot read textproto parameters %s", FLAGS.proto)
      return False

  # Overrides from the command line.
  # If the proto does not have a value, just use the default from FLAGS.
  if not options.proto.HasField("eta"):
    options.proto.eta = FLAGS.eta

  if not options.proto.HasField("max_depth"):
    options.proto.max_depth = FLAGS.max_depth

  if not options.proto.HasField("n_estimators"):
    options.proto.n_estimators = FLAGS.n_estimators

  if not options.proto.HasField("subsample"):
    options.proto.subsample = FLAGS.subsample

  if not options.proto.HasField("colsample_bytree"):
    options.proto.colsample_bytree = FLAGS.colsample_bytree

  if not options.proto.HasField("colsample_bylevel"):
    options.proto.colsample_bylevel = FLAGS.colsample_bylevel

  if not options.proto.HasField("colsample_bynode"):
    options.proto.colsample_bynode = FLAGS.colsample_bynode
  if not options.proto.HasField("reg_alpha"):
    options.proto.reg_alpha = FLAGS.reg_alpha
  if not options.proto.HasField("reg_lambda"):
    options.proto.reg_lambda = FLAGS.reg_lambda
  if not options.proto.HasField("gamma"):
    options.proto.gamma = FLAGS.gamma


  if option_present("tree_method"):
    match FLAGS.tree_method:
      case "auto":
        options.proto.tree_method = xgboost_model_pb2.AUTO
      case "exact":
        options.proto.tree_method = xgboost_model_pb2.EXACT
      case "approx":
        options.proto.tree_method = xgboost_model_pb2.APPROX
      case "hist":
        options.proto.tree_method = xgboost_model_pb2.HIST
      case _:   # Cannot happen because this is DEFINE_enum
        print(f"Unrecognised tree method {FLAGS.tree_method}", file=sys.stderr)
        return False
  elif not options.proto.HasField("tree_method"):
    options.proto.tree_method = xgboost_model_pb2.AUTO

  if not build_xgboost_model(argv[1], FLAGS.activity, options):
    logging.error("Model %s not build", options.mdir)
    return False

  if FLAGS.rescore:
    rescore_training_set(options)

  # zero return code for success.
  return 0

if __name__ == '__main__':
  app.run(main)
