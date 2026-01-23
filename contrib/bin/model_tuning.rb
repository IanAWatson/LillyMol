#!/usr/bin/env ruby
# Build svmfp models based on a pre-split set of files
# and a set of fingerprints.

require_relative 'lib/iwcmdline.rb'

def usage
  $stderr << "Build svmfp models\n"
  $stderr << " -trpct <pct>         percent of dataset to use for training\n"
  $stderr << " -fp <fname>          file of fingerprints to try\n"
  $stderr << " -nsplit <nsplit>     number of splits to create\n"
  exit(1)
end

cl = IWCmdline.new('-v-trpct=ipos-nsplit=ipos-niter=ipos-A=sfile-fp=sfile-catboost=close-xgboost=close-xgboost_config=sfile-lightgbm=close-lightgbm_config=sfile-i=ipos-keep_models')

if cl.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage
end

verbose = cl.option_present('v')

if ARGV.empty?
  $stderr << "Must specify smiles file\n"
  usage
end

smiles = ARGV[0]

if ! cl.option_present('trpct')
  $stderr << "Must specify the training set percent via the -trpct option\n"
  usage
end

trpct = 80
prefix = "A#{trpct}"
if cl.option_present('trpct')
  x = cl.value('trpct')
  if x < 0 || x > 100
    $stderr << "The training percent (-trpct) option must be a valid percent\n"
    usage
  end
  prefix = "A#{x}"
end

lightgbm = cl.values('lightgbm')
lightgbm_config = cl.value('lightgbm_config')
catboost = cl.values('catboost')
xgboost = cl.values('xgboost')
xgboost_config = cl.value('xgboost_config')

if ! cl.option_present('fp')
  $stderr << "Must specify file of fingerprints via the -fp option\n"
  usage
end

fingerprints = File.readlines(cl.value('fp'))

if ! cl.option_present('A')
  $stderr << "Must specify activity file via the -A option\n"
  usage
end

activity_fname = cl.value('A')

def execute_cmd(cmd, verbose, expected_files)
  $stderr << "Executing '#{cmd}'\n" if verbose
  system(cmd)
  expected_files.each do |fname|
    if ! File.size?(fname)
      "$stderr << "#{cmd} did not create #{fname}\n"
      return false
    end
  end
  return true
end

# The -ps option has been given. We need to determine the files implied by this.
# Return a list of the train and test splits.
def get_pre_split(stem)
  ndx = 0
  loop do
    train_smi = "TRAIN#{ndx}.smi"
    test_smi = "TEST#{ndx}.smi"
    return ndx if File.size?(train_smi) && File.size?(test_smi)
    ndx += 1
    return 0 if ndx > 1000
  end
end

if cl.option_present('ps')
  train_files, test_files = get_pre_split(cl.value('ps'))
  nsplit = train_files.size
else
  if cl.option_present('nsplit')
    nsplit = cl.value('nsplit')
  elif cl.option_present('niter')
    nsplit = cl.value('niter')
  else
    nsplit = 10
  end

  cmd = "stratified_samples -s 1 -N #{nsplit} -p #{trpct} -R TRAIN -E TEST -M #{smiles} #{activity_fname}"
  execute_cmd(cmd, verbose, ['TRAIN0.smi', 'TEST0.smi'])
  train_files = (0..nsplit).map { |i| "TRAIN#{i}.smi"}
  test_files = (0..nsplit).map { |i| "TEST#{i}.smi"}
end

stem = if cl.option_present('stem')
    cl.value('stem')
  else
    'model'
  end

support = if cl.option_present('p')
    cl.value('p')
  else
    1
  end

keep_models = cl.option_present('keep_models')

cmd_stream = File.open("model_tuning.txt", "w")

iwstats = "iwstats -w -Y allequals -p 2 -E #{activity_fname}"

svmfp_make = "svmfp_make.sh -A #{activity_fname}"
svmfp_make << " -p #{support}" if support > 1

fingerprints.each do |fp|
  fp.chomp!
  fps = fp.gsub(/ /, '')
  preds = []
  (0...nsplit).each do |split|
    mdir = "#{stem}.#{fps}.#{split}"
    cmd_stream << "#{svmfp_make} -mdir #{mdir} -gfp #{fp} -gfp #{train_files[split]}\n"
    pred = "#{stem}_#{fps}_#{split}.pred"
    cmd_stream << "svmfp_evaluate.sh -mdir #{mdir} #{test_files[split]} > #{pred}\n"
    preds << pred
    cmd_stream << "#{iwstats} #{pred} > #{prefix}.#{fps}.#{split}\n"
    cmd_stream << "rm -r #{mdir}\n" unless keep_models
    catboost.each_with_index do |cb, ndx|
      c_mdir = "#{mdir}_catboost_#{ndx}"
      cmd_stream << "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -catboost #{cb} -catboost -A #{activity_fname} #{train_files[split]}\n"
      pred = "#{stem}_#{fps}_cb_#{split}_#{ndx}.pred"
      cmd_stream << "catboost_evaluate.sh -mdir #{c_mdir} #{test_files[split]} > #{pred}\n"
      preds << pred
      cmd_stream << "#{iwstats} #{pred} > #{prefix}.ctbfp.#{ndx}.#{split}\n"
      cmd_stream << "rm -r #{c_mdir}\n" unless keep_models
    end
    xgboost.each_with_index do |xg, ndx|
      c_mdir = "#{mdir}_xgboost"
      cmd_stream << "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -xgboost -xgboost -xgboost_config #{xgboost_config} -A #{activity_fname} #{train_files[split]}\n"
      pred = "#{stem}_#{fps}_xg_#{split}_#{ndx}.pred"
      cmd_stream << "xgboost_evaluate.sh -v -mdir #{c_mdir} #{test_files[split]} > #{pred}\n"
      preds << pred
      cmd_stream << "#{iwstats} #{pred} > #{prefix}.xgbfp.#{ndx}.#{split}\n"
      cmd_stream << "rm -r #{c_mdir}\n" unless keep_models
    end
    lightgbm.each_with_index do |lg, ndx|
      c_mdir = "#{mdir}_lightgbm_#{ndx}"
      cmd_stream << "svmfp_make.sh -mdir #{c_mdir} -gfp #{fp} -gfp -lightgbm #{lg} -lightgbm -lightgbm_config #{lightgbm_config} -A #{activity_fname} #{train_files[split]}\n"
      pred = "#{stem}_#{fps}_lg_#{split}_#{ndx}.pred"
      cmd_stream << "lightgbm_evaluate.sh -mdir #{c_mdir} #{test_files[split]} > #{pred}\n"
      preds << pred
      cmd_stream << "#{iwstats} #{pred} > #{prefix}.lgbfp.#{ndx}.#{split}\n"
      cmd_stream << "rm -r #{c_mdir}\n" unless keep_models
    end
  end
# ifile = "I#{stem}_#{fps}"
# cmd = "iwstats -E #{activity_fname} -p 2 -I #{ifile} #{preds.join(' ')}"
# execute_cmd(cmd, verbose, [ifile])
end
